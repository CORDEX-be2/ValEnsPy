# helper functions for unit conversion - can be moved to more appropriate place

## imports for helpen functions
import xarray as xr
import pandas as pd
import numpy as np
from valenspy._utilities import load_yml

CORDEX_VARIABLES = load_yml("CORDEX_variables")

def convert_all_units_to_CF(ds: xr.Dataset, raw_LOOKUP, metadata_info: dict):
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
    unit_conversion_functions = {
        "Celcius": _convert_Celcius_to_Kelvin,
        "hPa": _convert_hPa_to_Pa,
        "mm": _convert_mm_to_kg_m2s,
    }

    for raw_var in ds.data_vars:
        var = next(
            (k for k, v in raw_LOOKUP.items() if v.get("raw_name") == raw_var), None
        )

        if var:  # Dont processes variables that are not in the lookup table.
            
            # convert units - based on the raw units
            raw_units = raw_LOOKUP[var]["raw_units"]
            if raw_units in unit_conversion_functions:
                ds = ds.rename_vars({raw_var: var}) # rename variable to CORDEX variable name

                #Do the conversion
                ds[var] = unit_conversion_functions[raw_units](ds[var])

                #Metadata attributes
                ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var]["standard_name"]
                ds[var].attrs["long_name"] = CORDEX_VARIABLES[var]["long_name"]
                ds[var].attrs["original_name"] = raw_LOOKUP[var]["raw_name"]
                ds[var].attrs["original_long_name"] = raw_LOOKUP[var]["raw_long_name"]
                ds[var].attrs["original_units"] = raw_units
            
            ds[var]["time"] = pd.to_datetime(ds[var].time)

            if metadata_info:
                for key, value in metadata_info.items():
                    ds[var].attrs[key] = value
            
            #Be carefull here because there should be some type of warning that the units are not converted
            #This could be because they are in the the correct unit or because the unit conversion function is not implemented
            #In the latter case an incorrect renaming of the variable could have happened
                

    return ds

# Do we want other possible inputs than data arrays?
def _convert_Celcius_to_Kelvin(da: xr.DataArray):
    """
    Convert values in xarray DataArray from °C to K

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # do conversion
    da = da + 273.15  # Celcius to Kelvin

    # update units attribute --  naming of units defined in ./src/valenspy/ancilliary_data/CORDEX_variables.yml
    da.attrs["units"] = "K"

    return da


def _convert_Kelvin_to_Celcius(da: xr.DataArray):
    """
    Convert values in xarray DataArray from K to °C

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # do conversion
    da = da - 273.15  # Kelvin to Celcius

    # update units attribute
    da.attrs["units"] = "°C"

    return da


def _convert_hPa_to_Pa(da: xr.DataArray):
    """
    Convert values in xarray DataArray from hPa to Pa

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # do conversion
    da = da * 100

    # update units attribute
    da.attrs["units"] = "Pa"

    return da


def _convert_Pa_to_hPa(da: xr.DataArray):
    """
    Convert values in xarray DataArray from Pa to hPa

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # do conversion
    da = da / 100

    # update units attribute
    da.attrs["units"] = "hPa"

    return da


def _convert_mm_to_kg_m2s(da: xr.DataArray):
    """
    Convert daily (!) values in xarray DataArray from mm to kg m^-2 s^-1

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # first, get timestep (frequency) by calculating the difference between the first consecutive time values in seconds
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    # do conversion
    da = da / timestep_nseconds  # mm to kg m^-2 s^-1

    # update units attribute
    da.attrs["units"] = "kg m-2 s-1"

    return da


def _convert_m_to_kg_m2s(da: xr.DataArray):
    """
    Convert values in xarray DataArray from mm hr^-1 to kg m^-2 s^-1

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # do conversion
    da = da * 1000 / 3600  # mm hr^-1 to kg m^-2 s^-1

    # update units attribute
    da.attrs["units"] = "kg m-2 s-1"

    return da

    # Convert J/m²/hr to W/m²
    w_m2 = j_m2_hr / seconds_per_hour


def _convert_J_m2_to_W_m2(da: xr.DataArray):
    """
    Convert values in xarray DataArray from J m^2 to W m^2

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """
    # first, get timestep (frequency) by calculating the difference between the first consecutive time values in seconds
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    # do conversion
    da = da / timestep_nseconds  # J m^2 to W m^2

    # update units attribute
    da.attrs["units"] = "W m-2"

    return da


def _convert_kWh_m2_day_to_W_m2(da: xr.DataArray):
    """
    Convert values in xarray DataArray from kWh/m2/day to W m^2

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """

    # do conversion
    da = da * (1000) / 24

    # update units attribute
    da.attrs["units"] = "W m-2"

    return da