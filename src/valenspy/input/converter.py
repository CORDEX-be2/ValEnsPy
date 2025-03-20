from pathlib import Path
from typing import Callable, Union
import xarray as xr
from valenspy._utilities import load_xarray_from_data_sources, load_yml, cf_status
from valenspy.input.unit_converter import convert_all_units_to_CF
from valenspy.input.converter_functions import (
    EOBS_to_CF,
    ERA5_to_CF,
    CCLM_to_CF,
    ALARO_K_to_CF,
    RADCLIM_to_CF,
)


class InputConverter:
    """A class that converts input files to netCDF files in CF convention."""

    def __init__(self, lookup, converter: Callable = None, metadata_info: dict = None):
        """Initialize the InputProcessor with a converter function.

        Parameters
        ----------
        lookup : str | dict | file
            The lookup table to convert variables (names and units) to CF convention.
            It can be a string (name of the lookup table in the ancilliary_data folder),
            It can be a dictionary (keys are the variable names and values are the CF names),
            It can be a json file (path to the json file containing the lookup table).
        converter : function, optional
            A function that deals with unique aspects of the input data when converting to CF convention. Default is None.
            This function is applied after the units are converted to CF convention and is the last step in the conversion process.
        metadata_info : dict, optional
            A dictionary containing the metadata information for the netCDF file.
        """
        self.converter = converter
        self.lookup = lookup if isinstance(lookup, dict) else load_yml(lookup)
        self.metadata_info = metadata_info

    def convert_input(self, data_sources, metadata_info=None):
        """Convert the input file(s)/xarray dataset to CF convention.

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

        ds = convert_all_units_to_CF(ds, self.lookup, metadata_info)
        ds = _set_global_attributes(ds, metadata_info)
        
        cf_status(ds)

        return ds

def _set_global_attributes(ds: xr.Dataset, metadata_info):
    for key, value in metadata_info.items():
        ds.attrs[key] = value
    return ds

INPUT_CONVERTORS = {
    "ERA5": InputConverter("ERA5_lookup", ERA5_to_CF, metadata_info={"dataset": "ERA5"}),
    "ERA5-Land": InputConverter("ERA5_lookup", ERA5_to_CF, metadata_info={"dataset": "ERA5-Land"}),
    "EOBS": InputConverter("EOBS_lookup", EOBS_to_CF, metadata_info={"freq": "day", "spatial_resolution": "0.1deg", "region": "Europe", "dataset": "EOBS"}),
    "CLIMATE_GRID": InputConverter("CLIMATE_GRID_lookup", metadata_info={"freq": "day", "spatial_resolution": "0.07° x 0.045° (~5km)", "region": "Belgium", "dataset": "CLIMATE_GRID"}),
    "CCLM": InputConverter("CCLM_lookup", CCLM_to_CF, metadata_info={"dataset": "CCLM"}),
    "ALARO_K": InputConverter("ALARO-SFX_K_lookup", ALARO_K_to_CF, metadata_info={"dataset": "ALARO_K"}),
    "RADCLIM": InputConverter("RADCLIM_lookup", RADCLIM_to_CF, metadata_info={"freq": "hour", "region": "Belgium", "dataset": "RADCLIM"}),
}

# Idea is to extend the shared functionality here (with subclasses if required) while the inputconvertor_functions are model specific.

# Needed:
#  - Some helper functions to extend input to glob arguments, str arguments etc.
#  - Extend input so that already loaded datasets can also be input
#  - CF Checker functionality
