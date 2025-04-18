"""Defines the InputManager class for loading and managing input data for ValEnsPy."""
from pathlib import Path
import pandas as pd
import xarray as xr
from intake_esm.cat import Assets, Attribute, AggregationControl, ESMCatalogModel
import intake
from xarray import DataTree
import re
import os
import glob
import warnings


from valenspy.input.converter import INPUT_CONVERTORS
from valenspy._utilities import load_yml, create_named_regex

DATASET_PATHS = load_yml("dataset_PATHS")
CORDEX_VARIABLES = load_yml("CORDEX_variables")

#TODO - update documentation once re-implemented including type hints
#This is an type of ECGtools light version

class InputManager:
    COLS = [
    #Required for unique identification of the dataset
    "source_id", #The source_id is the name of the dataset (e.g. "ERA5", "CNRM-CM6-1")
    "activity_id", #The activity_id is the name of the activity (e.g. "CORDEX", "CMIP6", "reanalysis", "etc")
    "domain_id", #The domain_id is the name of the domain (e.g. "europe", "global", "etc")
    "frequency", #The frequency is the frequency of the data (e.g. "hourly", "daily", "monthly", "yearly")
    "resolution", #The resolution is the resolution of the data (e.g. "0.11", "0.44", "etc")
    "version", #The version is the version of the dataset (e.g. "v1", "v2", "etc")
    "experiment_id", #The experiment_id is the name of the experiment (e.g. "historical", "rcp85", "etc")
    "driving_source_id", #The driving_source_id is the name of the driving dataset (e.g. "ERA5", "CNRM-CM6-1") including variant label CNRM-CM6-1_r1i1p1f1
    "realization_id", #The realization_id is the realization of the dataset (e.g. "r1i1p1f1", "r2i1p1f1", "etc" or numbering)
    #File specific metadata
    "variable_id",
    "time_range",
    ]

    def __init__(self, machine : str, dataset_info : dict = None):
        """
        Initialize the InputManager.

        Parameters
        ----------
        machine : str
            The name of the machine (used to identify dataset paths).
        dataset_info : dict
            A dictionary containing dataset information. The keys are dataset names and the values are dictionaries with the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
                <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
            - meta_data: A dictionary containing metadata for the dataset.
            Default is None. If None, the built-in dataset_info for the provided machine is used. See the dataset_PATHS.yaml file.
        """
        self.machine = machine

        if dataset_info:
            self.datasets_yaml = dataset_info
        else:
            self.datasets_yaml = DATASET_PATHS[machine]

        self._validate_dataset_info()

        self.df = self.create_catalog()

    @classmethod
    def from_yaml(cls, yaml_path):
        """
        Create an InputManager instance from a ESM-catalog YAML file.
        """
        self.cat = intake.open_esm_datastore(yaml_path)
        self.df = self.cat._df

    @property
    def readable_catalog(self):
        """
        Return a readable version of the catalog DataFrame.
        """
        COLS = [
        #Required for unique identification of the dataset
        "source_id", #The source_id is the name of the dataset (e.g. "ERA5", "CNRM-CM6-1")
        "frequency", #The frequency of the data (e.g. "daily", "monthly")

        "variable_id",
        "time_range",
        ]
        cols = [item for item in COLS if item not in ["variable_id", "time_range"]]
        def count(x):
            return len(x)

        return self.df.groupby(cols)["variable_id"].apply(count).to_frame()

    @property
    def available_datasets(self):
        """
        Return a list of available datasets in the catalog.
        """
        return self.df["source_id"].unique().tolist()

    def create_intake_esm_json_catalog(
        self, 
        output_path,
        path_column_name="path",
        variable_column_name="variable_id",
        data_format="netcdf",
        groupby_attrs=None,
        aggregations=None,
        esmcat_version:str = "0.1.0",
        description:str=None,
    ):
        """
        Create an intake-esm JSON catalog from the datasets in the catalog.

        Inspired on ecgtools save functionality
        """
        for col in {variable_column_name, path_column_name}.union(set(groupby_attrs or [])):
            assert col in self.df.columns, f"Column {col} not found in DataFrame"

        attributes = [Attribute(column_name=column, vocabulary="") for column in self.df.columns]

        _aggregation_control = AggregationControl(
            variable_column_name=variable_column_name,
            groupby_attrs=groupby_attrs,
            aggregations=aggregations,
        )

        cat = ESMCatalogModel(
            esmcat_version=esmcat_version,
            description=description,
            attributes=attributes,
            aggregation_control=_aggregation_control,
            assets=Assets(column_name=path_column_name, format=data_format),
        )

        cat._df = self.df

        self.cat = cat

    def _validate_dataset_info(self):

        required_identifiers = [
            "source_id",
            "activity_id",
            "domain_id"
        ]

        required_for_filtering = [
            "frequency",
            "period"
        ]

        for dataset_name, dataset_info in self.datasets_yaml.items():
            key_set = set(dataset_info.get("meta_data", {}).keys())
            key_set.update({"source_id"}) 
            pattern = dataset_info.get("pattern", None)
            if pattern:
                # Extract the identifiers from the pattern
                key_set.update(re.findall(r"<(.*?)>", pattern))

            # Check if all required identifiers are present
            for identifier in required_identifiers:
                if identifier not in key_set:
                    warnings.warn(f"Dataset {dataset_name} is missing the required identifier '{identifier}' in its pattern or metadata.")

            # Check if all required identifiers for filtering are present
            for identifier in required_for_filtering:
                if identifier not in key_set:
                    warnings.warn(f"Dataset {dataset_name} is missing the required identifier '{identifier}' for filtering in its pattern or metadata.")

    def update_catalog_from_yaml(self, yaml_path):
        """
        Update the catalog from a YAML file.

        Parameters
        ----------
        yaml_path : Path
            The path to the YAML file containing dataset information.
        """
        # Load the YAML file
        datasets_info = load_yml(yaml_path)
        for dataset_name, dataset_info in datasets_info.items():
            # Add the new dataset to the catalog
            self.update_catalog(dataset_name, dataset_info)
        
    def update_catalog_from_dataset_info(self, dataset_name, dataset_root_dir, dataset_pattern, metadata={}):
        """
        Update the catalog with a new dataset.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset.
        dataset_root_dir : str
            The root directory of the dataset.
        dataset_pattern : str
            The regex pattern for matching files in the dataset.
        metadata : dict, optional
            Additional metadata to include in the catalog (default is empty dictionary).
        """
        dataset_info = {
            "root": dataset_root_dir,
            "pattern": dataset_pattern,
            "meta_data": metadata,
        }
        self.update_catalog(dataset_name, dataset_info)
    
    def update_catalog(self, dataset_name, dataset_info_dict):
        """
        Add a new dataset to the catalog.

        Parameters
        ----------
        dataset_info_dict : dict
            A dictionary containing dataset information. The keys are dataset names and the values are dictionaries with the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset.
            - meta_data: A dictionary containing metadata for the dataset.
        """
        self.datasets_yaml[dataset_name] = dataset_info_dict
        self._validate_dataset_info()
        data = self._process_dataset_for_catalog(dataset_name, self.datasets_yaml[dataset_name])
        df = pd.DataFrame(data)
        self.df = pd.concat([self.df, df], ignore_index=True)

    def _process_dataset_for_catalog(self, dataset_name, dataset_info):
        """
        Process all files in a dataset and extract metadata.
        """
        dataset_root = Path(dataset_info.get("root"))
        regex_pattern = create_named_regex(dataset_info.get("pattern", None))
        regex = re.compile(dataset_root.as_posix() + r"/" + regex_pattern)

        dataset_meta_data = dataset_info.get("meta_data", {})

        IC = INPUT_CONVERTORS.get(dataset_name, None)
        if IC:
            CORDEX_variable_set = IC.cordex_variables
            variable_set = IC.raw_variables
            long_name_set = IC.raw_variables_long_names

        files_with_metadata = []
        for root, _, files in os.walk(dataset_root):
            for file in files:
                if file.endswith(".nc"):
                    file_path = os.path.join(root, file)
                    if match := regex.match(file_path):
                        file_metadata = match.groupdict()
                    else:
                        file_metadata = {}

                    # Add the file path to the metadata
                    file_metadata["path"] = Path(file_path)

                    # Add dataset level metadata
                    file_metadata = {**dataset_meta_data, **file_metadata}
                    file_metadata["source_id"] = dataset_name

                    # Translate the variable_id to the CORDEX variable name (if possible)
                    if IC:
                        variable_id = file_metadata.get("variable_id")
                        if not variable_id:
                            file_metadata["raw_variable_id"] = list(variable_set)
                            file_metadata["variable_id"] = list(CORDEX_variable_set)
                        elif variable_id in variable_set or variable_id in long_name_set:
                            file_metadata["raw_variable_id"] = variable_id
                            file_metadata["variable_id"] = IC.get_CORDEX_variable(variable_id)

                    #Convert time data to a time range
                    if "year" in file_metadata:
                        start_year = file_metadata["year"]
                        end_year = file_metadata["year"]
                    elif "start_year" in file_metadata and "end_year" in file_metadata:
                        start_year = file_metadata["start_year"]
                        end_year = file_metadata["end_year"]
                    elif "yearmonthday" in file_metadata:
                        start_year = file_metadata["yearmonthday"][:4]
                        end_year = file_metadata["yearmonthday"][:4]
                    else:
                        start_year = dataset_info.get("start_year", None)
                        end_year = dataset_info.get("end_year", None)
                    if start_year and end_year:
                        file_metadata["start_year"] = start_year
                        file_metadata["end_year"] = end_year
                        try:
                            file_metadata["time_range"] = pd.Interval(
                                left=pd.Timestamp(f"{start_year}-01-01"),
                                right=pd.Timestamp(f"{end_year}-12-31"),
                                closed="both"
                            )
                        except Exception as e:
                            file_metadata["time_range"] = None
                    else:
                        file_metadata["start_year"] = None
                        file_metadata["end_year"] = None
                        file_metadata["time_range"] = None
                    
                    files_with_metadata.append(file_metadata)

        return files_with_metadata
        
    def create_catalog(self):
        """
        Create a catalog by scanning dataset paths and extracting metadata.
        """
        files_with_metadata = []
        for dataset_name, dataset_info in self.datasets_yaml.items():
            # Process the dataset and extract metadata
            grouped_files_with_metadata = self._process_dataset_for_catalog(dataset_name, dataset_info)
            # Add the dataset name to the metadata
            files_with_metadata.extend(grouped_files_with_metadata)
            
        # Create a DataFrame and save it as a CSV
        df = pd.DataFrame(files_with_metadata)
        return df

    def open_dataset(
        self,
        dataset_name,
        variables=["tas"],
        period=None,
        freq=None,
        other_filters={},
        cf_convert=True,
        metadata_info={},
    ):
        """
        Load a dataset and return an xarray DataArray or Dataset.
        """

        df = self.df

        number_of_files = {"catalog": len(df)}

        # Check if the dataset name is valid
        if dataset_name not in self.available_datasets:
            raise ValueError(f"Dataset {dataset_name} not found in catalog. Available datasets: {self.available_datasets}")

        df = df[df["source_id"] == dataset_name]
        # Filter the DataFrame based on the provided parameters
        if variables:
            def filter_variables(x):
                if isinstance(x, str):
                    return x in variables
                elif isinstance(x, list):
                    return any(var in variables for var in x)
                else:
                    return False
            df = df[df["variable_id"].apply(filter_variables)]
            number_of_files["variables"] = len(df)
        if period:
            pass
        if freq:
            df = df[df["frequency"] == freq]
            number_of_files["frequency"] = len(df)
        if other_filters:
            for key, value in other_filters.items():
                if key in df.columns:
                    df = df[df[key] == value]
                    number_of_files[key] = len(df)
        
        if df.empty:
            raise ValueError(f"No data found for dataset {dataset_name} with the specified filters.\n Filters: {variables}, {period}, {freq}, {other_filters} \n Number of files per filter: {number_of_files}")
        
        # Check if an input converter is available for the dataset
        IC = INPUT_CONVERTORS.get(dataset_name, None)
        if IC and cf_convert:
            ds = IC.convert_input(df["path"].to_list(), metadata_info=metadata_info)
        else:
            ds = xr.open_mfdataset(df["path"], decode_coords="all", chunks="auto")

        return ds

    def open_datatree(
        self,
        dataset_paths,
        variables=["tas"],
        period=None,
        freq=None,
        other_filters={},
        cf_convert=True,
        metadata_info={},
    ):
        """
        Load multiple datasets and return a DataTree of xarray DataArrays or Datasets.

        Note the dataset_paths 
        """
        datatree_dict = {}
        
        for dataset_path in dataset_paths:
            dataset_name = dataset_path.split("/")[-1]
            if dataset_name not in self.available_datasets:
                raise ValueError(f"Dataset {dataset_name} not found in catalog. Available datasets: {self.available_datasets}")
            else:
                datatree_dict[dataset_path] = self.open_dataset(
                    dataset_name,
                    variables=variables,
                    period=period,
                    freq=freq,
                    other_filters=other_filters,
                    cf_convert=cf_convert,
                    metadata_info=metadata_info,
                )
        return DataTree.from_dict(datatree_dict)

