"""Defines the InputConverter class for converting input files to ValEnsPy complaint xarrays."""

from pathlib import Path
from typing import Callable, Union
import xarray as xr
from valenspy._utilities import (
    load_xarray_from_data_sources,
    load_yml,
    cf_status,
    _set_global_attributes,
    _convert_all_units_to_CF
)
from valenspy.input.converter_functions import (
    EOBS_to_CF,
    ERA5_to_CF,
    CCLM_to_CF,
    ALARO_K_to_CF,
    RADCLIM_to_CF,
    MAR_to_CF
)

class InputConverter:
    """A class for converting input files or xarrays to ValEnsPy complaint xarrays."""

    def __init__(self, var_lookup_table: str | dict | Path, converter: Callable = None, metadata_info: dict = None):
        """Initialize the InputProcessor

        Parameters
        ----------
        var_lookup_table : str | dict | Path
            A dictionary, or a Path to a yml file or a string matching the name of one of the yml files in valenspy/ancillary_data, 
            Keys are the CORDEX standard variable names and values are information about the variable in the input file.
        converter : function, optional
            A function that deals with unique aspects of the input data when converting to CF convention. Default is None.
            This function is applied before the units and variable names are converted.
        metadata_info : dict, optional
            A dictionary containing the metadata information for the netCDF file. This is added to the global attributes of the netCDF file.

        Examples
        --------
        >>> from valenspy.input.converter import InputConverter
        >>> ERA5_dict = {
        ...     "tas": {
        ...         "raw_name": "t2m",
        ...         "raw_units": "K"}
        ...  }
        >>> converter = InputConverter(var_lookup_table=ERA5_dict)
        >>> ds = converter([paths_to_era5_files])

        Or use the pre-defined input converters

        >>> from valenspy.input.converter import InputConverter
        >>> converter = InputConverter("ERA5_lookup")
        >>> ds = converter([paths_to_era5_files])
        """
        self.converter = converter
        self.var_lookup_table = var_lookup_table if isinstance(var_lookup_table, dict) else load_yml(var_lookup_table)
        self.metadata_info = metadata_info

    def __call__(self, data_sources: Path | list[Path] | xr.Dataset, metadata_info: dict = {}) -> xr.Dataset:
        """Convert the input file(s) or xarray dataset to CF convention."""
        return self.convert_input(data_sources, metadata_info)

    def convert_input(self, data_sources: Path | list[Path] | xr.Dataset, metadata_info: dict = {}) -> xr.Dataset:
        """Convert the input file(s) or xarray dataset to CF convention.

        Parameters
        ----------
        data_sources : Path or list(Path) or xarray.Dataset
            The input file or list of input files or an xarray dataset to convert.
        metadata_info : dict, optional
            A dictionary containing additional metadata information for the netCDF file.

        Returns
        -------
        xarray.Dataset
            An xarray dataset in CF convention.
        """
        ds = load_xarray_from_data_sources(data_sources)
        
        if self.converter:
            ds = self.converter(ds)
        
        metadata_info = {**self.metadata_info, **metadata_info}

        ds = _convert_all_units_to_CF(ds, self.var_lookup_table, metadata_info)
        ds = _set_global_attributes(ds, metadata_info)
        
        cf_status(ds)

        return ds

INPUT_CONVERTORS = {
    "ERA5": InputConverter("ERA5_lookup", 
                           ERA5_to_CF, 
                           metadata_info={"dataset": "ERA5"}),
    "ERA5-Land": InputConverter("ERA5_lookup", 
                                ERA5_to_CF, 
                                metadata_info={"dataset": "ERA5-Land"}),
    "EOBS": InputConverter("EOBS_lookup", 
                           EOBS_to_CF, 
                           metadata_info={"freq": "day", "spatial_resolution": "0.1deg", "region": "Europe", "dataset": "EOBS"}),
    "CLIMATE_GRID": InputConverter("CLIMATE_GRID_lookup", 
                                   metadata_info={"freq": "day", "spatial_resolution": "0.07° x 0.045° (~5km)", "region": "Belgium", "dataset": "CLIMATE_GRID"}),
    "CCLM": InputConverter("CCLM_lookup", 
                           CCLM_to_CF, 
                           metadata_info={"dataset": "CCLM"}),
    "ALARO_K": InputConverter("ALARO-SFX_K_lookup", 
                              ALARO_K_to_CF,
                              metadata_info={"dataset": "ALARO_K"}),
    "RADCLIM": InputConverter("RADCLIM_lookup", 
                              RADCLIM_to_CF, 
                              metadata_info={"freq": "hour", "region": "Belgium", "dataset": "RADCLIM"}),
    "MAR": InputConverter("MAR_lookup", 
                          MAR_to_CF, 
                          metadata_info={"dataset": "MAR", "freq": "day", "region": "Belgium"}),
}