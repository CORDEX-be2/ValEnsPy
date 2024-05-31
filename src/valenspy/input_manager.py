#Load the dataset_PATHS.yml
from yaml import safe_load
from pathlib import Path
from itertools import permutations
import xarray as xr

from valenspy.inputconverter import INPUT_CONVERTORS

src_path = Path(__file__).resolve().parent

with open(src_path / "ancilliary_data" / "dataset_PATHS.yml") as file:
    DATASET_PATHS = safe_load(file)

with open(src_path / "ancilliary_data" / "CORDEX_variables.yml") as file:
    CORDEX_VARIABLES = safe_load(file)

class InputManager:
    def __init__(self, machine):
        self.machine = machine
        self.dataset_paths = DATASET_PATHS[machine]

    def load_data(self, dataset_name, variables=["tas"], period=None, freq="daily", cf_convert=True, path_identifiers=[]):
        """
        Load the data for the specified dataset, variables, period and frequency and transform it into ValEnsPy CF-Compliant format.
        A regex is used to match any file paths that start with the dataset_path and contain the variables, period and frequency.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to load. This should be in the dataset_PATHS.yml file for the specified machine.
        variables : list, optional
            The variables to load. The default is ["tas"]. These should be CORDEX variables defined in CORDEX_variables.yml.
        period : list, optional
            The period to load start and end dates. The default is None.
        freq : str, optional
            The frequency of the data. The default is "daily".
        cf_convert : bool, optional
            Whether to convert the data to CF-Compliant format. The default is True.
        path_identifiers : list, optional
            Other identifiers to match in the file paths. These are on top the variable long name, year and frequency. The default is [].

        Returns
        -------
        ds : xarray.Dataset
            The loaded dataset in CF-Compliant format.
        """
        if self._is_valid_dataset_name(dataset_name):
            files = self._get_file_paths(dataset_name, variables=variables, period=period, freq=freq, path_identifiers=path_identifiers)
            if cf_convert:
                input_converter = INPUT_CONVERTORS[dataset_name]
                ds = input_converter.convert_input(files)
            else:
                ds = xr.open_mfdataset(files)
        return ds
           

    def _get_file_paths(self, dataset_name, variables=["tas"], period=None, freq="daily", path_identifiers=[]):
        """Get the file paths for the specified dataset, variables, period and frequency."""
        with open(src_path / "ancilliary_data" / f"{dataset_name}_lookup.yml") as file:
            obs_LOOKUP = safe_load(file)
        dataset_path = Path(self.dataset_paths[dataset_name])
        file_paths = []
        for variable in variables:
            obs_long_name = obs_LOOKUP[variable]["obs_long_name"]
            if period:
                for year in range(period[0], period[1]+1):
                    for pattern in generate_regex([obs_long_name, year, freq] + path_identifiers):
                        file_paths+=dataset_path.glob(pattern)
        return set(file_paths)

    def _is_valid_dataset_name(self, dataset_name):
        """Check if the dataset name is valid for the machine."""
        if not dataset_name in self.dataset_paths:
            raise ValueError(f"Dataset name {dataset_name} is not valid for machine {self.machine}. Valid dataset names are {list(self.dataset_paths.keys())}. See dataset_PATHS.yml.")
        return True

def generate_regex(components):
    perm = permutations([str(c) for c in components])    
    # Create a regex pattern that matches any of the permutations containing the components and ending in .nc
    patterns = [f"**/*{'*'.join(p)}*.nc" for p in perm]
    return patterns