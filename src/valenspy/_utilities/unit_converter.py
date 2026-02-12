# helper functions for unit conversion - can be moved to more appropriate place

## imports for helpen functions
import xarray as xr
import numpy as np
import xclim
import pandas as pd
from valenspy._utilities import load_yml

CORDEX_VARIABLES = load_yml("CORDEX_variables")

def _convert_all_units_to_CF(ds: xr.Dataset, raw_LOOKUP, metadata_info: dict):
    """Convert all units for all variables in the dataset to the correct units by applying the correct conversion function.

    The raw_units attribute is used to determine the current units and therefore the conversion function to apply.
    The target units are defined in the CORDEX_variables.yml file. In addition it:
    - Handels variable metadata attributes specific to the var and conversion
    - Adds dataset metadata attributes to the variable attributes
    - Converts the time dimension to a pandas datetime object


    Parameters
    ----------
    ds : xr.Dataset
        The xarray dataset to convert
    raw_LOOKUP : dict
        The lookup table for the variables which corresponds the CORDEX variables to the variables in the dataset and their units
    metadata_info : dict
        The metadata information which should be added to the each variable in the dataset

    Returns
    -------
    xr.Dataset
        The converted xarray dataset
    """

    for var,var_lookup in raw_LOOKUP.items():

        raw_var = var_lookup.get("raw_name")

        if raw_var in ds:
            raw_units = var_lookup.get("raw_units")

            standard_name = CORDEX_VARIABLES[var].get("standard_name")
            long_name = CORDEX_VARIABLES[var].get("long_name")
            units = CORDEX_VARIABLES[var].get("units")

            #Rename the variable
            ds = ds.rename_vars({raw_var: var})

            #Add units 
            if "units" not in ds[var].attrs:
                ds[var].attrs["units"] = raw_units
            #Correct units if they are not the same as the raw units
            elif ds[var].attrs["units"] != raw_units:
                ds[var].attrs["ds_original_units"] = ds[var].attrs["units"]
                ds[var].attrs["units"] = raw_units

            #Add the standard name for help with CF unit conversions see https://xclim.readthedocs.io/en/stable/notebooks/units.html#Smart-conversions:-Precipitation
            if "standard_name" not in ds[var].attrs and "raw_standard_name" in var_lookup:
                ds[var].attrs["standard_name"] = var_lookup.get("raw_standard_name")
            elif "standard_name" in ds[var].attrs:
                if "raw_standard_name" in var_lookup:
                    ds[var].attrs["ds_original_standard_name"] = ds[var].attrs["standard_name"]
                    ds[var].attrs["standard_name"] = var_lookup.get("raw_standard_name")
            
            #Use xclim to handle all unit conversions from the raw units (ds[var].attrs["units"]) to the target units 
            #units attribute is automatically updated
            ds[var] = xclim.units.convert_units_to(ds[var], units, context="infer")

            #Overwrite the xclim standard name with the CORDEX standard name
            ds[var].attrs["standard_name"] = standard_name

            # Add metadata
            ds[var].attrs["long_name"] = long_name
            ds[var].attrs["original_name"] = raw_var
            ds[var].attrs["original_units"] = raw_units

            ds[var]["time"] = pd.to_datetime(ds[var].time)

            if metadata_info:
                for key, value in metadata_info.items():
                    ds[var].attrs[key] = value

            if "freq" not in ds[var].attrs:
                ds[var].attrs["freq"] = xr.infer_freq(ds[var].time)
            
    return ds