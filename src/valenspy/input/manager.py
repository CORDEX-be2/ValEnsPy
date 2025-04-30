"""Defines the InputManager class for finding, managing, preprocessing and loading input data for ValEnsPy."""
from functools import partial
from pathlib import Path

from valenspy.input.converter import INPUT_CONVERTORS
from valenspy.input.catalog import ValenspyEsmDatastore
from valenspy.input.esm_catalog_builder import CatalogBuilder, CATALOG_COLS
from valenspy._utilities import load_yml

esmcat_default_data = {
    "id" : "test",
    "esmcat_version": "0.1.0",  # intake-esm JSON file structure version, as per: https://github.com/NCAR/esm-collection-spec
    "assets": {"column_name": "path", "format": "netcdf"},
    "aggregation_control": {
        "variable_column_name": "variable_id",
        "groupby_attrs": CATALOG_COLS["required_identifiers"] + CATALOG_COLS["required_identifiers_with_default"],
        "aggregations": [
            {
                "type": "join_existing",
                "attribute_name": "time_period_start",
                "options": {"dim": "time"},
            },
            {
                "type": "union",   
                "attribute_name": "variable_id"},
        ],
    },
    "attributes": [],
}

#TODO: Inhert from both ValenspyEsmDatastore and CatalogBuilder? - this may allow:
#  - default values for to_dataset_dict methods for preprocess, xarray_open_kwargs and xarray_combine_by_coords_kwargs
#  - No need for the __getattr__ method
class InputManager:
    """A class to find, manage, preprocess and load input data for ValEnsPy.
    
    The InputManager class consists of an ValEnsPy specific intake-esm catalog (ValenspyEsmDatastore) and a CatalogBuilder. 
    The Catalog Builder is used to create the catalog, a df with dataset information per file, using minimal information about the datasets and their path structure. 
    This catalog is then used to create an esm_datastore (ValenspyEsmDatastore) which can be used to search and load the datasets.
    The InputManager class provides a preprocessing function based on the input convertors to convert the datasets to ValEnsPy  
    """

    def __init__(
            self, 
            machine : str, 
            datasets_info : dict = None, 
            description=None, 
            input_convertors : dict = INPUT_CONVERTORS,
            esmcat_data : dict = esmcat_default_data,
            xarray_open_kwargs : dict = {},
            xarray_combine_by_coords_kwargs : dict = {},
            intake_esm_kwargs : dict = {"sep": "/"}
            ):
        """
        Initialize an InputManager.

        Parameters
        ----------
        machine : str
            The name of the catalog. If dataset_info is not passed it will be used to load the dataset_info from the built-in dataset_paths.yaml file.
        datasets_info : dict | str, optional
            A dictionary containing datasets and their dataset information needed to build the catalog. This can be a dictionary or a path to a YAML file.
            Default is None. If None, the built-in dataset info for the provided catalog_id is used if it exists.
            The dictionary should contain dataset names as keys and their dataset information as values.
            The datasetinfo should contain the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
                <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
            - meta_data: A dictionary containing metadata for the dataset.
        description : str, optional
            A description of the catalog. Default is None. This is used to create the description in the intake-esm catalog.
        input_convertors : dict, optional
            A dictionary containing input convertors for the datasets. The keys are dataset names and the values are functions that convert the dataset to ValEnsPy compliant format.
            Default is INPUT_CONVERTORS, the set of predefined input convertors.
        esmcat_data : dict, optional
            A dictionary containing the description of esmcat data to create the intake-esm catalog. Default is esmcat_default_data. See :func:`intake_esm.esm_datastore.from_dict` for more information.
        xarray_open_kwargs : dict, optional
            A dictionary containing default arguments to pass to xarray_open_kwargs in :func:`intake_esm.esm_datastore.to_dataset_dict`, :func:`intake_esm.esm_datastore.to_dask` and/or :func:`intake_esm.esm_datastore.to_datatree`.
        xarray_combine_by_coords_kwargs : dict, optional
            A dictionary containing default arguments to pass to xarray_combine_by_coords_kwargs in :func:`intake_esm.esm_datastore.to_dataset_dict`, :func:`intake_esm.esm_datastore.to_dask` and/or :func:`intake_esm.esm_datastore.to_datatree`.
        intake_esm_kwargs : dict, optional
            A dictionary containing additional arguments for the creation of the intake_esm catalog. Default is an empty dictionary. See :func:`intake_esm.esm_datastore` for more information.
        """
        self.catalog_builder = CatalogBuilder(
            catalog_id=machine,
            datasets_info=datasets_info
        )

        if self.catalog_builder.skipped_files:
            print("There are files that were skipped during the catalog creation. Please check the skipped_files property for a list of skipped files.")

        #Collection of functions to preprocess the data so that they are ValEnsPy compliant
        self.input_convertors = input_convertors

        if description:
            esmcat_data["description"] = description

        self.esm_datastore = ValenspyEsmDatastore(
            obj={"esmcat": esmcat_data, "df": self.catalog_builder.df},
            **intake_esm_kwargs
            )

        #Options for xarray.open_dataset and xarray.combine_by_coords
        self.xarray_open_kwargs = xarray_open_kwargs
        self.xarray_combine_by_coords_kwargs = xarray_combine_by_coords_kwargs

    def __getattr__(self, name):
        """
        Delegate attribute access to the esm_datastore instance (self.esm_datastore).

        This allows all methods and attributes of esm_datastore to be accessed
        directly from the InputManager instance.
        """
        if hasattr(self.esm_datastore, name):
            return getattr(self.esm_datastore, name)
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    # @classmethod
    # def from_yaml(cls, yaml_path, intake_esm_kwargs : dict = {}):
    #     """
    #     Create an InputManager instance from a YAML file.

    #     Parameters
    #     ----------
    #     yaml_path : str
    #         The path to the YAML file containing dataset information. It should contain:
    #         - root: The root directory of the dataset.
    #         - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
    #             <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
    #         - meta_data: A dictionary containing metadata for the dataset.
    #     intake_esm_kwargs : dict
    #         A dictionary containing additional arguments for the intake_esm catalog. Default is an empty dictionary. See the intake_esm documentation for more information.
    #     """
    #     # Load the YAML file
    #     datasets_info = load_yml(yaml_path)
    #     # Create an instance of InputManager
    #     return cls(machine=None, dataset_info=datasets_info, intake_esm_kwargs=intake_esm_kwargs)
    
    @property
    def skipped_files(self):
        """
        The files that where skipped during the catalog creation.
        """
        return self.catalog_builder.skipped_files
    
    @property
    def intake_to_xarray_kwargs(self):
        """
        Easy access of kwargs to be used passed that can be passed to :func:`intake_esm.esm_datastore.to_dataset_dict`, :func:`intake_esm.esm_datastore.to_dask` and/or :func:`intake_esm.esm_datastore.to_datatree`

        Three types of kwargs are created:
        - preprocess: The preprocessor function which applies the input convertor to the dataset if an input convertor exists (i.e. source_id is in INPUT_CONVERTORS).
        - xarray_open_kwargs: The kwargs to be passed to xarray.open_dataset.
        - xarray_combine_by_coords_kwargs: The kwargs to be passed to xarray.combine_by_coords.
        """
        kwargs = {
            "preprocess": self.preprocess,}
        if self.xarray_open_kwargs:
            kwargs["xarray_open_kwargs"] = self.xarray_open_kwargs
        if self.xarray_combine_by_coords_kwargs:
            kwargs["xarray_combine_by_coords_kwargs"] = self.xarray_combine_by_coords_kwargs
        return kwargs

    def add_input_convertor(self, dataset_name, input_convertor):
        """Add an input convertor to the InputManager."""
        self.input_convertors[dataset_name] = input_convertor

    def _update_catalog(self, dataset_name, dataset_info_dict):
        """
        Update the catalog (df with dataset information per file) with a new dataset.

        The catalog_builder is used to parse the dataset information and update the catalog.

        Parameters
        ----------
        dataset_info_dict : dict
            A dictionary containing dataset information. The keys are dataset names and the values are dictionaries with the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset.
            - meta_data: A dictionary containing metadata for the dataset.
        """
        self.catalog_builder.add_dataset(
            dataset_name,
            dataset_info_dict
        )

    def update_catalog_from_yaml(self, yaml_path):
        """
        Update the catalog from a YAML file.

        For each dataset, parse the dataset information, validate it, add it to the catalog and update the esm_datastore.

        Parameters
        ----------
        yaml_path : Path
            The path to the YAML file containing datasets information a dictionary of dataset names and their dataset information.
            The datasetinfo should contain the following keys:
            - root: The root directory of the dataset.
            - pattern: The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
                <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
            - meta_data: A dictionary containing metadata for the dataset.
        """
        # Load the YAML file
        datasets_info = load_yml(yaml_path)
        for dataset_name, dataset_info in datasets_info.items():
            # Add the new dataset to the catalog
            self._update_catalog(dataset_name, dataset_info)

        self.catalog_builder._validate_dataset_info()
        self.esm_datastore.esmcat._df = self.catalog_builder.df
        
    def update_catalog_from_dataset_info(self, dataset_name, dataset_root_dir, dataset_pattern, metadata={}):
        """
        Update the catalog with a new dataset.

        For the dataset, parse the dataset information, validate it, add it to the catalog and update the esm_datastore.

        Parameters
        ----------
        dataset_name : str
            The name of the dataset.
        dataset_root_dir : str
            The root directory of the dataset.
        dataset_pattern : str
            The regex pattern for matching files in the dataset. This is the reletave path starting from the root and in the following format:
            <indentifier_name>/<indentifier_name>/<indentifier_name>_fixed_part_<variable_id>/<another_identifier>_<year>.nc
        metadata : dict, optional
            Additional metadata to include in the catalog. Default is an empty dictionary.
        """
        dataset_info = {
            "root": dataset_root_dir,
            "pattern": dataset_pattern,
            "meta_data": metadata,
        }
        self._update_catalog(dataset_name, dataset_info)
        self.catalog_builder._validate_dataset_info()
        self.esm_datastore.esmcat._df = self.catalog_builder.df

    @property
    def preprocess(self):
        """
        A preprocessor function to convert the input dataset to ValEnsPy compliant data.

        This function applys the input convertor to the dataset if an input convertor exists (i.e. source_id is in this managers input convertors).
        """

        def process_IC(ds, IC_dict, df):
            file_name = ds.encoding["source"]
            source_id = df[df["path"] == Path(file_name)]["source_id"].values[0]
            if source_id in IC_dict:
                return IC_dict[source_id](ds)
            else:
                return ds
            
        return partial(process_IC, IC_dict=self.input_convertors, df=self.esm_datastore.df)