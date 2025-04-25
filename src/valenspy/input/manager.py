"""Defines the InputManager class for loading and managing input data for ValEnsPy."""
from pathlib import Path
import pandas as pd
import xarray as xr
from intake_esm import esm_datastore
from intake_esm.cat import Assets, Attribute, Aggregation, AggregationControl, ESMCatalogModel
from xarray import DataTree
import re
import os
import warnings
import copy


from valenspy.input.converter import INPUT_CONVERTORS
from valenspy._utilities import load_yml, create_named_regex, parse_time_period

DATASET_PATHS = load_yml("dataset_PATHS")
CORDEX_VARIABLES = load_yml("CORDEX_variables")

#TODO - update documentation once re-implemented including type hints
#This is an type of ECGtools light version

CATALOG_COLS = {
#Required for unique identification of the dataset
# - each unique combination of these columns should uniquely define an xarray dataset
"required_identifiers" :[
    "source_id", #The source_id is the name of the dataset (e.g. "ERA5", "CNRM-CM6-1")
    "source_type", #The source_type is the type of the dataset (e.g. "reanalysis", "CMIP6", "CMPI5-CORDEX", "observations")
    "domain_id", #The domain_id is the name of the domain (e.g. "europe", "EUR-11", "global", "etc")
    "experiment_id", #The experiment_id is the name of the experiment (e.g. "historical", "rcp85", "ssp585", "land_use_change", "etc")
    "version", #The version is the version of the dataset (e.g. "v1", "v2", "etc")
    "resolution", #The resolution is the resolution of the data (e.g. "0.11", "0.44", "etc")
    "frequency", #The frequency is the frequency of the data (e.g. "hourly", "daily", "monthly", "yearly")
],
#Required but default values are used if not relevant
"required_identifiers_with_default" : [
    "driving_source_id", #The driving_source_id is the name of the driving dataset (e.g. "ERA5", "CNRM-CM6-1") including variant label CNRM-CM6-1_r1i1p1f1
    "institution_id", #The institution_id is the name of the institution (e.g. "CNRM", "KMI", "KNMI", "etc")
    "realization", #The realization_id is the realization of the dataset (e.g. "r1i1p1f1", "r2i1p1f1", "etc" or numbering)
],
"filtering_identifiers" : [
    #Filtering identifiers within a unique dataset allowing to limit the number of files to load
    "variable_id",
    "time_period_start", #Note that time_period will be created from time_period or time_period_start/time_period_end and time_format
    "time_period_end",
]
}

#TODO: Check if this should not be a direct subclass of esm_datastore
class InputManager:
    """A class to manage input data for ValEnsPy using intake-esm."""

    def __init__(self, machine : str, dataset_info : dict = None, input_convertors : dict = INPUT_CONVERTORS, intake_esm_kwargs : dict = {}):
        """
        Initialize an InputManager.

        Parameters
        ----------
        machine : str
            The name of the machine. If dataset_info is not passed it will be used to load the dataset_info from the built-in dataset_paths.yaml file.
        dataset_info : dict
            A dictionary containing dataset information. The keys are dataset names and the values are dictionaries with the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
                <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
            - meta_data: A dictionary containing metadata for the dataset.
            Default is None. If None, the built-in dataset_info for the provided machine is used. See the dataset_PATHS.yaml file.
        intake_esm_kwargs : dict
            A dictionary containing additional arguments for the intake_esm catalog. Default is an empty dictionary. See the intake_esm documentation for more information.
        """
        self.machine = machine
        self.input_convertors = input_convertors

        if dataset_info:
            self.datasets_yaml = dataset_info
        else:
            self.datasets_yaml = DATASET_PATHS[machine]

        self._validate_dataset_info()

        self.skipped_files = {}

        self.df = self.create_df()

        self.intake_kwargs = intake_esm_kwargs

        self.create_intake_esm_json_from_df(**self.intake_kwargs)

    #All other functions applied on the manager should be applied on the catalog
    def __getattr__(self, name):
        """
        Delegate attribute access to the esm_datastore instance (self.esm_datastore).

        This allows all methods and attributes of esm_datastore to be accessed
        directly from the InputManager instance.
        """
        if hasattr(self.esm_datastore, name):
            return getattr(self.esm_datastore, name)
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    @classmethod
    def from_yaml(cls, yaml_path, intake_esm_kwargs : dict = {}):
        """
        Create an InputManager instance from a YAML file.

        Parameters
        ----------
        yaml_path : str
            The path to the YAML file containing dataset information. It should contain:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
                <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
            - meta_data: A dictionary containing metadata for the dataset.
        intake_esm_kwargs : dict
            A dictionary containing additional arguments for the intake_esm catalog. Default is an empty dictionary. See the intake_esm documentation for more information.
        """
        # Load the YAML file
        datasets_info = load_yml(yaml_path)
        # Create an instance of InputManager
        return cls(machine=None, dataset_info=datasets_info, intake_esm_kwargs=intake_esm_kwargs)

    def add_input_convertor(self, dataset_name, input_convertor):
        """
        Add an input convertor to the InputManager.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset.
        input_convertor : object
            An instance of the input convertor class.
        """
        self.input_convertors[dataset_name] = input_convertor

    @property
    def available_datasets(self):
        """
        Return a list of available datasets in the catalog.
        """
        return self.df["source_id"].unique().tolist()

    def create_intake_esm_json_from_df(
        self, 
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

        Parameters
        ----------
        path_column_name : str
            The name of the column containing the file paths. Default is "path".
        variable_column_name : str
            The name of the column containing the variable names. Default is "variable_id".
        data_format : str
            The format of the data files. Default is "netcdf".
        groupby_attrs : list
            A list of attributes to group the data by. The 
        aggregations : list
            A list of intake-esm Aggregation attributes to aggregate the data. The default (None) aggregates as follows:
            - "union" on the variable_column_name
            - "join_existing" on the time_period attribute

        Inspired on ecgtools save functionality
        """
        for col in {variable_column_name, path_column_name}.union(set(groupby_attrs or [])):
            assert col in self.df.columns, f"Column {col} not found in DataFrame"

        #Possibly something is wrong with this aggregation type causing issues for the esm_datastore when trying to call to_dataset_dict
        if aggregations is None:
            aggregations = [
                Aggregation(type="union", 
                            attribute_name=variable_column_name), 
                Aggregation(type="join_existing", 
                            attribute_name="time_period", 
                            options={
                                "dim": "time",
                                "coords": "minimal",
                                "compat": "override"
                            }
                    ),
                ]

        if groupby_attrs is None:
            groupby_attrs = CATALOG_COLS["required_identifiers"] + CATALOG_COLS["required_identifiers_with_default"]

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

        self.esm_datastore = esm_datastore(cat)

    def _validate_dataset_info(self):

        required_identifiers = CATALOG_COLS["required_identifiers"]
        required_identifiers_with_default = CATALOG_COLS["required_identifiers_with_default"]
        filtering_identifiers = CATALOG_COLS["filtering_identifiers"]

        for dataset_name, dataset_info in self.datasets_yaml.items():
            key_set = set(dataset_info.get("meta_data", {}).keys())
            pattern = dataset_info.get("pattern", None)
            if pattern:
                for key in re.findall(r"<(.*?)>", pattern):
                    key_set.add(key)

            # Check if all required identifiers are present
            for identifier in required_identifiers:
                if identifier not in key_set:
                    warnings.warn(f"Dataset {dataset_name} is missing the required identifier '{identifier}' in its pattern or metadata.")

            # Check if all required identifiers with default values are present
            for identifier in required_identifiers_with_default:
                if identifier not in key_set:
                    # Set default value if not present
                    dataset_meta_data = dataset_info.get("meta_data", {})
                    dataset_meta_data[identifier] = "default_value"

            # Check if all required identifiers for filtering are present
            for identifier in filtering_identifiers:
                if identifier == "time_period_start" or identifier == "time_period_end": #This is checked seperately
                    continue
                if identifier not in key_set:
                    warnings.warn(f"Dataset {dataset_name} is missing the required identifier '{identifier}' for filtering in its pattern or metadata.")

            # Check if time_period or time_period_start/time_period_end is present
            if ("time_period" not in key_set) and ("time_period_start" not in key_set and "time_period_end" not in key_set):
                warnings.warn(f"Dataset {dataset_name} is missing the required identifier 'time_period' or 'time_period_start/time_period_end' in its pattern or metadata.")

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
        data = self._process_dataset_for_catalog(dataset_name, self.datasets_yaml[dataset_name])
        df = pd.DataFrame(data)
        self.df = pd.concat([self.df, df], ignore_index=True)

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
        self._validate_dataset_info()
        self.create_intake_esm_json_from_df()
        
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
        self._validate_dataset_info()
        self.create_intake_esm_json_from_df()

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
                        #Add the skipped file to the skipped files dictionary (create the entry if it does not exist)
                        if dataset_name not in self.skipped_files:
                            self.skipped_files[dataset_name] = []
                        self.skipped_files[dataset_name].append(file_path)
                        continue

                    # Add the file path to the metadata
                    file_metadata["path"] = Path(file_path)

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

                    #Create time_period attribute using time_period or time_period_start/time_period_end and the time_format
                    time_period = file_metadata.pop("time_period", None)
                    time_period_start = file_metadata.pop("time_period_start", None)
                    time_period_end = file_metadata.pop("time_period_end", None)
                    time_format = file_metadata.pop("time_format", None)                  

                    if time_period and not (time_period_start or time_period_end):
                        start, end = parse_time_period(time_period, format=time_format)
                    elif time_period_start and time_period_end:
                        start,_ = parse_time_period(time_period_start, format=time_format) #The earliest timestamp of that period, e.g. 2024 -> 2024-01-01 00:00:00
                        _,end = parse_time_period(time_period_end, format=time_format) #The latest timestamp of that period, e.g. 2024 -> 2024-12-31 23:59:59
                    else:
                        start, end = None, None

                    if start and end:
                        file_metadata["time_period_start"] = start
                        file_metadata["time_period_end"] = end
                    else:
                        self.skipped_files[dataset_name].append(file_path)
                        continue

                    files_with_metadata.append(file_metadata)  

        return files_with_metadata
        
    def create_df(self):
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

    #Some hybrid esm-intake functionality as esm-intake to_dataset_dict does not seem to work with the current catalog for mulitple files
    def open_dataset(
        self,
        preprocessor = True,
        open_xarray_kwargs = {"decode_coords":"all", "chunks":"auto"},
        **query
    ):
        """
        Load a single dataset from the catalog and return an xarray DataArray or Dataset.

        Parameters
        ----------
        preprocessor : bool
            If True, apply the input converter to the dataset if an input converter exists (i.e. source_id is in INPUT_CONVERTORS). Default is True.
        **query : 
            keyword arguments to filter the catalog. keywords should be the same as the columns in the catalog.
        """

        subcat = self.search(**query)
        if len(subcat.keys()) > 1:
            raise ValueError(f"Multiple datasets found for the given query. Please specify a more specific query. Found datasets: {subcat.keys()}")

        elif len(subcat.keys()) == 0:
            raise ValueError(f"No datasets found for the given query. Explore the catalog with self.search() to find available datasets.")

        else:
            paths_column = subcat.esmcat.assets.column_name
            files = subcat.df[paths_column].to_list()

            dataset_name = subcat.df["source_id"].unique()[0] #There is only one unique source_id as we only have one dataset

            metadata_info = subcat.keys_info().to_dict('records')[0]

            return self._open_xarray_dataset(
                dataset_name,
                files,
                metadata_info,
                preprocessor=preprocessor,
                **open_xarray_kwargs
            )
        
    def search(self, require_all_on=None, **query):
        esm_datastore = copy.deepcopy(self.esm_datastore)
        
        if "time_period" in query.keys():
            df = self.esmcat.df.copy()
            if isinstance(query["time_period"], str):
                start, end = parse_time_period(query["time_period"])
            elif isinstance(query["time_period"], list):
                start, _ = parse_time_period(query["time_period"][0])
                _, end = parse_time_period(query["time_period"][1])
            else:
                raise ValueError("time_period should be a string or a list of strings")

            #Filter keeping files which cover a period which overlaps with the time_period
            df = df[(pd.to_datetime(df["time_period_start"]) <= end) & (start <= pd.to_datetime(df["time_period_end"]))]
            esm_datastore.esmcat._df = df
            query.pop("time_period")

        return esm_datastore.search(
            require_all_on=require_all_on,
            **query
        )

    def open_datatree(
        self,
        preprocessor = True,
        source_id_extra_queries = None,
        tree_structure = None,
        open_xarray_kwargs = {"decode_coords":"all", "chunks":"auto"},
        **query
    ):
        """
        Create a DataTree from a search query on the catalog. Each node is a unique dataset defined by the groupby_attrs in the catalog which has relevant data based on the query.
        """

        subcat = self.search(**query)

        datatree_dict = {}

        for key in subcat.keys():
            dataset_name = subcat[key].df["source_id"].unique()[0]
            metadata_info = subcat.keys_info().to_dict('records')[0]
            files = subcat[key].df["path"].to_list()
            if not tree_structure:
                path = key.replace(".", "/")
            else:
                if tree_structure[0] == "/":
                    tree_structure = tree_structure[1:]
                    path = "/"
                else:
                    path = ""
                path += "/".join(metadata_info[path] for path in tree_structure.split("/"))
            
            datatree_dict[path] = self._open_xarray_dataset(
                dataset_name,
                files,
                metadata_info,
                preprocessor=preprocessor,
                **open_xarray_kwargs
            )
        
        return DataTree.from_dict(datatree_dict)

    def _open_xarray_dataset(
            self,
            dataset_name,
            files,
            metadata_info,
            preprocessor = True,
            **kwargs
    ):
        """
        Open an xarray dataset from a list of files and metadata information.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset.
        files : list
            A list of file paths to open.
        metadata_info : dict
            A dictionary containing metadata information for the dataset.
        preprocessor : bool
            If True, apply the input converter to the dataset if an input converter exists (i.e. source_id is in INPUT_CONVERTORS). Default is True.
        """
        # Check if an input converter is available for the dataset
        IC = self.input_convertors.get(dataset_name, None)
        if IC and preprocessor:
            ds = IC(files, metadata_info=metadata_info)
        else:
            ds = xr.open_mfdataset(files, **kwargs)
        return ds