"""Defines the InputManager class for loading and managing input data for ValEnsPy."""
from pathlib import Path
import xarray as xr
from datatree import DataTree
import re
import glob

from valenspy.input.converter import INPUT_CONVERTORS
from valenspy._utilities import load_yml

DATASET_PATHS = load_yml("dataset_PATHS")
CORDEX_VARIABLES = load_yml("CORDEX_variables")

#TODO - update documentation once re-implemented including type hints

class InputManager:
    def __init__(self, machine):
        self.machine = machine
        self.dataset_paths = DATASET_PATHS[machine]

    def load_m_data(
        self, datasets_dict, variables=["tas"], cf_convert=True, metadata_info={}
    ):
        """
        Load multiple datasets and variables and return a DataTree object.

        Each dataset is passed to the load_data method and the resulting datasets are combined into a DataTree object.

        Parameters
        ----------
        datasets_dict : dict
            A dictionary of datasets to load. The keys are the dataset names and the values are dictionaries containing the period, frequency,
            region and path_identifiers as keys.
        variables : list
            The variables to load. The default is ["tas"]. These should be CORDEX variables defined in CORDEX_variables.yml.
        cf_convert : bool, optional
            Whether to convert the data to CF-Compliant format. The default is True.
        metadata_info : dict, optional
            Other metadata information to pass to the input converter. The default is {}.

        Returns
        -------
        DataTree
            A DataTree object containing the loaded datasets.

        Examples
        --------
        >>> manager = InputManager(machine='hortense')
        >>> # Get all ERA5 tas (temperature at 2m) at a daily frequency for the years 2000 and 2001. The paths must include "max".
        >>> data_request_dict={
            "EOBS":
                {"path_identifiers":["mean"]},
            "ERA5":
                {"period":[2000,2001],
                "freq":"daily",
                "region":"europe",
                "path_identifiers":["min"]}
            }
        >>> dt = manager.load_m_data(data_request_dict, variables=["tas","pr"])
        """

        ds_dict = {}
        for dataset_name, dataset_info in datasets_dict.items():
            # pass all the dataset info as kwargs to the load_data method
            print(f"Loading data for {dataset_name}...")
            ds_dict[dataset_name] = self.load_data(
                dataset_name,
                variables=variables,
                cf_convert=cf_convert,
                metadata_info=metadata_info,
                **dataset_info,
            )
        return DataTree.from_dict(ds_dict)

    def load_data(
        self,
        dataset_name,
        variables=["tas"],
        period=None,
        freq=None,
        region=None,
        cf_convert=True,
        path_identifiers=[],
        metadata_info={},
    ):
        """
        Load the data for the specified dataset, variables, period and frequency and transform it into ValEnsPy CF-Compliant format.

        For files to be found and loaded they should be in a subdirectory of the dataset path and contain
        the raw_long_name or raw_name or CORDEX variable name, the year (optional), frequency and path_identifiers (optional) in the file name.

        A regex search is used to match any netcdf (.nc) file paths that start with the dataset_path from the dataset_PATHS.yml and contains:
        1) The raw_long_name of the CORDEX variables given the dataset_name_lookup.yml
        2) Any YYYY string within the period
        3) The frequency of the data (daily, monthly, yearly)
        4) Any additional path_identifiers

        The order of these components is irrelevant. The dataset is then loaded using xarray.open_mfdataset and if cf_convert is True, the data is converted
        to CF-Compliant format using the appropriate input converter. If no period is specified, all files matching the other components are loaded.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to load. This should be in the dataset_PATHS.yml file for the specified machine.
        variables : list, optional
            The variables to load. The default is ["tas"]. These should be CORDEX variables defined in CORDEX_variables.yml.
        period : list or an int, optional
            The period to load. If a list, the start and end years of the period. For a single year both an int and a list with one element are valid. The default is None.
        freq : str, optional
            The frequency of the data. The default is None.
        region : str, optional
            The region to load. The default is None.
        cf_convert : bool, optional
            Whether to convert the data to CF-Compliant format. The default is True.
        path_identifiers : list, optional
            Other identifiers to match in the file paths. These are on top the variable long name, year and frequency. The default is [].
        other_metadata_info : dict, optional
            Other metadata information to pass to the input converter. The default is {}.

        Returns
        -------
        ds : xarray.Dataset
            The loaded dataset in CF-Compliant format.

        Raises
        ------
        FileNotFoundError
            If no files are found for the specified dataset, variables, period, frequency and path_identifiers.

        ValueError
            If the dataset name is not valid for the machine. i.e. not in the dataset_PATHS.yml file.

        Examples
        --------
        >>> manager = InputManager(machine='hortense')
        >>> # Get all ERA5 tas (temperature at 2m) at a daily frequency for the years 2000 and 2001. The paths must include "max".
        >>> ds = manager.load_data("ERA5", variables=["tas"], period=[2000,2001], path_identifiers=["max"])
        """
        if isinstance(period, list):
            if len(period) > 2:
                raise ValueError("Period must be a list at most 2 elements or an int.")
            if len(period) == 1:
                period = int(period[0])

        if self._is_valid_dataset_name(dataset_name):
            files = self._get_file_paths(
                dataset_name,
                variables=variables,
                period=period,
                freq=freq,
                region=region,
                path_identifiers=path_identifiers,
            )
            if not files:
                raise FileNotFoundError(
                    f"No files found for dataset {dataset_name}, variables {variables}, period {period}, frequency {freq}, region {region} and path_identifiers {path_identifiers}."
                )
            print("File paths found:")
            for f in files:
                print(f)
            if cf_convert:
                input_converter = INPUT_CONVERTORS[dataset_name]
                if period:
                    metadata_info["period"] = period
                if freq:
                    metadata_info["freq"] = freq
                if region:
                    metadata_info["region"] = region
                if path_identifiers:
                    metadata_info["path_identifiers"] = path_identifiers
                ds = input_converter.convert_input(files, metadata_info=metadata_info)
            else:
                ds = xr.open_mfdataset(files, chunks="auto")
        return ds

    def _get_file_paths(
        self,
        dataset_name,
        variables=["tas"],
        period=None,
        freq=None,
        region=None,
        path_identifiers=[],
    ):
        """Get the file paths for the specified dataset, variables, period and frequency."""

        # ERA5Land has same lookuptable as ERA5
        if dataset_name == "ERA5-Land":
            dataset_name_lookup = "ERA5"
        else:
            dataset_name_lookup = dataset_name

        raw_LOOKUP = load_yml(f"{dataset_name_lookup}_lookup")

        dataset_path = Path(self.dataset_paths[dataset_name])
        file_paths = []
        variables = (
            [variables] if isinstance(variables, str) else variables
        )  # if single variable inputted as string, convert to list
        for variable in variables:
            if variable not in raw_LOOKUP:
                var_regex = f"{variable}"
            else:
                raw_long_name = raw_LOOKUP[variable]["raw_long_name"]
                raw_name = raw_LOOKUP[variable]["raw_name"]
                var_regex = f"({raw_long_name}|{raw_name}_|{variable}_)"
            components = [var_regex] + path_identifiers
            if period:
                if isinstance(period, int):
                    year_regex = f"({period})"
                else:
                    year_regex = f"({'|'.join([str(year) for year in range(period[0], period[1]+1)])})"
                components.append(year_regex)
            if freq:
                components.append(freq)
            if region:
                components.append(region)
            file_paths += [
                f
                for f in dataset_path.glob("**/*.nc")
                if all(
                    re.search(f"{dataset_path}/.*{component}.*", str(f))
                    for component in components
                )
            ]

        return list(set(file_paths))

    def _is_valid_dataset_name(self, dataset_name):
        """Check if the dataset name is valid for the machine."""
        if not dataset_name in self.dataset_paths:
            raise ValueError(
                f"Dataset name {dataset_name} is not valid for machine {self.machine}. Valid dataset names are {list(self.dataset_paths.keys())}. See dataset_PATHS.yml."
            )
        return True
