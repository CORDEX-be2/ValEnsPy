from pathlib import Path
from typing import Callable, Union
import xarray as xr
from valenspy._utilities import load_xarray_from_data_sources
from valenspy.input.converter_functions import (
    EOBS_to_CF,
    ERA5_to_CF,
    ERA5Land_to_CF,
    CLIMATE_GRID_to_CF,
    CCLM_to_CF,
    ALARO_K_to_CF,
    RADCLIM_to_CF,
)


class InputConverter:
    """A class that converts input files to netCDF files in CF convention."""

    def __init__(self, converter: Callable):
        """Initialize the InputProcessor with a converter function.

        Parameters
        ----------
        converter : function
            The function to convert the input file(s) to a netCDF file in CF convention.
            The function should take a string or list of strings as input and return a netCDF file.
        """
        self.converter = converter

    def convert_input(self, data_sources, metadata_info=None):
        """Convert the input file(s)/xarray dataset to CF convention.

        Parameters
        ----------
        data_sources : Path or list(Path) or xarray.Dataset
            The input file or list of input files or an xarray dataset to convert.

        Returns
        -------
        xarray.Dataset
            An xarray dataset in CF convention.
        """
        ds = load_xarray_from_data_sources(data_sources)
        return self.converter(ds, metadata_info)


INPUT_CONVERTORS = {
    "ERA5": InputConverter(ERA5_to_CF),
    "ERA5-Land": InputConverter(ERA5Land_to_CF),
    "EOBS": InputConverter(EOBS_to_CF),
    "CLIMATE_GRID": InputConverter(CLIMATE_GRID_to_CF),
    "CCLM": InputConverter(CCLM_to_CF),
    "ALARO_K": InputConverter(ALARO_K_to_CF),
    "RADCLIM": InputConverter(RADCLIM_to_CF),
}

# Idea is to extend the shared functionality here (with subclasses if required) while the inputconvertor_functions are model specific.

# Needed:
#  - Some helper functions to extend input to glob arguments, str arguments etc.
#  - Extend input so that already loaded datasets can also be input
#  - CF Checker functionality
