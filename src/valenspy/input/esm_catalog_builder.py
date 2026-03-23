from pathlib import Path
import pandas as pd
import re
import os
import warnings

from valenspy.input.converter import INPUT_CONVERTORS
from valenspy._utilities import load_yml, create_named_regex, parse_string_to_time_period

DATASET_PATHS = load_yml("dataset_info")
CORDEX_VARIABLES = load_yml("CORDEX_variables")

CATALOG_COLS = {
#Required for unique identification of the dataset
# - each unique combination of these columns should uniquely define an xarray dataset
"required_identifiers" :[
    "source_id",
    "source_type",
    "domain_id",
    "experiment_id",
    "version",
    "resolution",
    "frequency",
],
#Required but default values are used if not relevant
"required_identifiers_with_default" : [
    "driving_source_id",
    "institution_id",
    "realization",
    "post_processing",
],
"filtering_identifiers" : [
    #Filtering identifiers within a unique dataset allowing to limit the number of files to load
    "variable_id",
    "time_period_start", #Note that time_period_start and time_period_end will be created from time_period or (time_period_start and time_period_end) and time_format
    "time_period_end",
]
}

class CatalogBuilder:

    def __init__(self, catalog_id, datasets_info : dict | str = None, datasets_subset = None):
        """
        Initialize the CatalogBuilder with a catalog ID and dataset information.

        Parameters
        ----------
        catalog_id : str
            The ID of the catalog. If dataset_info is not provided, this will be used to load pre-existing dataset_info from ValEnsPy if it exists.
        datasets_info : dict | str, optional
            A dictionary containing datasets and their dataset information needed to build the catalog. This can be a dictionary or a path to a YAML file.
            Default is None. If None, the built-in dataset info for the provided catalog_id is used if it exists.
            The dictionary should contain dataset names as keys and their dataset information as values.
            The datasetinfo should contain the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
                <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
            - meta_data: A dictionary containing metadata for the dataset
        dataset_to_load : str | list, optional
            The name of the dataset(s) to load. If not provided, all datasets in the dataset_info will be loaded.
        """
        self.catalog_id = catalog_id

        if datasets_info:
            if isinstance(datasets_info, str):
                datasets_info = load_yml(datasets_info)
            self.datasets_info = datasets_info
        else:
            self.datasets_info = DATASET_PATHS[catalog_id]

        if datasets_subset:
            if isinstance(datasets_subset, str):
                datasets_subset = [datasets_subset]
            # Filter the datasets_info to only include the datasets to load
            self.datasets_info = {dataset_name: self.datasets_info[dataset_name] for dataset_name in datasets_subset if dataset_name in self.datasets_info}

        self._validate_dataset_info()

        self.skipped_files = {}

        self.df = self.create_df()

    def add_dataset(self, dataset_name, dataset_info):
        """
        Update the dataset information for a specific dataset.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to update.
        dataset_info : dict
            The new dataset information to update.
        """
        #TODO - check if the dataset_name is already in the datasets_info
        self.datasets_info[dataset_name] = dataset_info
        new_data = self._process_dataset_for_catalog(dataset_name, dataset_info)
        
        #Only continue if new data was found
        if new_data:
            # Cast the variable_id column as a list - Hack needed for intake_esm catalog to garantee has_multiple_variable_assets in its current version
            new_df = pd.DataFrame(new_data)
            new_df = new_df.explode(["variable_id", "raw_variable_id"], ignore_index=True)
        
            self.df = pd.concat([self.df, new_df], ignore_index=True)

    def _validate_dataset_info(self):
        """Validate the dataset information to ensure all required identifiers are present."""
        required_identifiers = CATALOG_COLS["required_identifiers"]
        filtering_identifiers = CATALOG_COLS["filtering_identifiers"]

        for dataset_name, dataset_info in self.datasets_info.items():
            key_set = set(dataset_info.get("meta_data", {}).keys())
            pattern = dataset_info.get("pattern", None)
            if pattern:
                for key in re.findall(r"<(.*?)>", pattern):
                    key_set.add(key)

            # Check if all required identifiers are present
            for identifier in required_identifiers:
                if identifier not in key_set:
                    warnings.warn(f"Dataset {dataset_name} is missing the required identifier '{identifier}' in its pattern or metadata.")

            # Check if all required identifiers for filtering are present
            for identifier in filtering_identifiers:
                if identifier == "time_period_start" or identifier == "time_period_end": #This is checked seperately
                    continue
                if identifier not in key_set:
                    warnings.warn(f"Dataset {dataset_name} is missing the required identifier '{identifier}' for filtering in its pattern or metadata.")

            # Check if time_period or time_period_start/time_period_end is present
            if ("time_period" not in key_set) and ("time_period_start" not in key_set and "time_period_end" not in key_set):
                warnings.warn(f"Dataset {dataset_name} is missing the required identifier 'time_period' or 'time_period_start/time_period_end' in its pattern or metadata.")

    def create_df(self):
        """
        Create a catalog by scanning dataset paths and extracting metadata.
        """
        files_with_metadata = []
        for dataset_name, dataset_info in self.datasets_info.items():
            # Process the dataset and extract metadata
            print(f"Processing dataset: {dataset_name}")
            grouped_files_with_metadata = self._process_dataset_for_catalog(dataset_name, dataset_info)
            # Add the dataset name to the metadata
            files_with_metadata.extend(grouped_files_with_metadata)
            
        # Create a DataFrame and save it as a CSV
        df = pd.DataFrame(files_with_metadata)

        # Explode the variable_id column as a list
        df = df.explode(["variable_id", "raw_variable_id"], ignore_index=True)

        return df
    
    def _process_dataset_for_catalog(self, dataset_name, dataset_info):
        """
        Process all files in a dataset and extract metadata for each file.

        Given a dataset name and its information, this function parses every file in the dataset, returning the metadata parsed from the file name
        and the dataset level metadata.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset to process.
        dataset_info : dict
            The dataset information containing the root directory, regex pattern, and metadata.

        Returns
        -------
        list
            A list of dictionaries containing metadata for each file in the dataset.

        """
        dataset_root = Path(dataset_info.get("root"))
        regex_pattern = create_named_regex(dataset_info.get("pattern", None))
        regex = re.compile(dataset_root.as_posix() + r"/" + regex_pattern)

        dataset_meta_data = dataset_info.get("meta_data", {})
        if dataset_IC := dataset_info.get("input_convertor", None):
            IC = INPUT_CONVERTORS.get(dataset_IC, None) #Use the specified ICs
        else:
            IC = INPUT_CONVERTORS.get(dataset_name, None) #Use the dataset name to find the IC
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
                        #Add the skipped file to the skipped files dictionary (create the entry if it does not exist)
                        if dataset_name not in self.skipped_files:
                            self.skipped_files[dataset_name] = []
                        self.skipped_files[dataset_name].append(file_path)
                        continue

                    # Add the file path to the metadata
                    file_metadata["path"] = Path(file_path).as_posix()

                    # Add dataset level metadata
                    file_metadata = {**dataset_meta_data, **file_metadata}

                    # Translate the variable_id to the CORDEX variable name (if possible)
                    if IC:
                        variable_id = file_metadata.get("variable_id")
                        if not variable_id:
                            file_metadata["raw_variable_id"] = list(variable_set)
                            file_metadata["variable_id"] = list(CORDEX_variable_set)
                        elif variable_id in variable_set or variable_id in long_name_set:
                            file_metadata["raw_variable_id"] = variable_id
                            file_metadata["variable_id"] = IC.get_CORDEX_variable(variable_id)
                    else:
                        #<variable_id> is already in the pattern and in CORDEX variable name convention
                        file_metadata["raw_variable_id"] = file_metadata.get("variable_id", None)

                    #Create time_period attribute using time_period or time_period_start/time_period_end and the time_format
                    time_period = file_metadata.pop("time_period", None)
                    time_period_start = file_metadata.pop("time_period_start", None)
                    time_period_end = file_metadata.pop("time_period_end", None)
                    time_format = file_metadata.pop("time_format", None)                  

                    if time_period and not (time_period_start or time_period_end):
                        start, end = parse_string_to_time_period(time_period, format=time_format)
                    elif time_period_start and time_period_end:
                        start,_ = parse_string_to_time_period(time_period_start, format=time_format) #The earliest timestamp of that period, e.g. 2024 -> 2024-01-01 00:00:00
                        _,end = parse_string_to_time_period(time_period_end, format=time_format) #The latest timestamp of that period, e.g. 2024 -> 2024-12-31 23:59:59
                    else:
                        start, end = None, None

                    if start and end:
                        file_metadata["time_period_start"] = start
                        file_metadata["time_period_end"] = end
                    else:
                        if dataset_name not in self.skipped_files:
                            self.skipped_files[dataset_name] = []
                        self.skipped_files[dataset_name].append(file_path)
                        continue

                    for key in CATALOG_COLS["required_identifiers_with_default"]:
                        if key not in file_metadata:
                            file_metadata[key] = "default"

                    files_with_metadata.append(file_metadata)

        if not files_with_metadata:
            warnings.warn(f"No valid files found for dataset {dataset_name}; \n Please check the dataset root {dataset_root} and pattern {dataset_info.get('pattern', None)}")

        return files_with_metadata
