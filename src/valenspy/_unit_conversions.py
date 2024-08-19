# helper functions for unit conversion - can be moved to more appropriate place

## imports for helpen functions
import xarray as xr
import pandas as pd
import numpy as np
import warnings
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
        "mm/hr": _convert_m_to_kg_m2s,
        "m": _convert_m_to_kg_m2s,
        "m/hr": _convert_m_to_kg_m2s,
        "J/m^2": _convert_J_m2_to_W_m2,
        "kWh/m2/day": _convert_kWh_m2_day_to_W_m2,
    }

    # Key: The unit of the raw data
    # Value: The unit of the CORDEX equivalent unit or the unit that is used to identify the conversion function
    EQUIVALENT_UNITS = {"degC": "Celcius", "m/s": "m s-1"}

    for raw_var in ds.data_vars:
        var = next(
            (k for k, v in raw_LOOKUP.items() if v.get("raw_name") == raw_var), None
        )

        if var:  # Dont processes variables that are not in the lookup table.

            # convert units - based on the raw units
            raw_units = raw_LOOKUP[var]["raw_units"]

            # If the raw_units are in the equivalent_units, use the replacement unit
            if raw_units in EQUIVALENT_UNITS:
                raw_units = EQUIVALENT_UNITS[raw_units]

            if raw_units in unit_conversion_functions:
                ds = ds.rename_vars(
                    {raw_var: var}
                )  # rename variable to CORDEX variable name
                ds[var] = unit_conversion_functions[raw_units](
                    ds[var]
                )  # Do the conversion

            elif raw_units == CORDEX_VARIABLES[var]["units"]:
                # If the raw_units are the same as the target units, just rename the variable
                ds = ds.rename_vars({raw_var: var})

            else:
                # Throw a warning that the raw_unit in the lookup table is not implemented
                cordex_var_units = CORDEX_VARIABLES[var]["units"]
                warnings.warn(
                    f"Unit conversion for {raw_units} to {cordex_var_units} is not implemented for variable {var}."
                )

            if var in ds:  # If the renaming occured add the metadata attributes
                # Metadata attributes
                ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var]["standard_name"]
                ds[var].attrs["long_name"] = CORDEX_VARIABLES[var]["long_name"]
                ds[var].attrs["units"] = CORDEX_VARIABLES[var]["units"]
                ds[var].attrs["original_name"] = raw_LOOKUP[var]["raw_name"]
                ds[var].attrs["original_long_name"] = raw_LOOKUP[var]["raw_long_name"]
                ds[var].attrs["original_units"] = raw_LOOKUP[var]["raw_units"]

                ds[var]["time"] = pd.to_datetime(ds[var].time)

                if metadata_info:
                    for key, value in metadata_info.items():
                        ds[var].attrs[key] = value

                # If freq is not in the metadata_info, we can try to infer it from the time dimension
                if "freq" not in ds[var].attrs:
                    ds[var].attrs["freq"] = _determine_time_interval(ds[var])

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


def _convert_kg_m2s_to_mh(da: xr.DataArray):
    """
    Convert values in xarray DataArray from kg m^-2 s^-1 to mm hr^-1

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
    da = da * 3600 / 1000  # kg m^-2 s^-1 to mm hr^-1

    # update units attribute
    da.attrs["units"] = "mm hr-1"

    return da


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


# helper functions - can be moved to more appropriate place


def _determine_time_interval(da: xr.DataArray):
    """
    Find the time interval (freq) of the input data array based on it's time axis, by calculating the difference between the first two time instances.

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray with time axis to check the time interval

    Returns
    -------
    freq : string
        The frequency string containing "hourly, daily, monthly or yearly"
    """

    diff = da.time.diff(dim="time").values[0]

    # Check for exact differences
    if diff == np.timedelta64(1, "h"):
        freq = "hourly"
    elif diff == np.timedelta64(1, "D"):
        freq = "daily"
    elif diff == np.timedelta64(1, "M"):
        freq = "monthly"
    elif diff == np.timedelta64(1, "Y"):
        freq = "yearly"
    else:
        return (
            "Difference does not match exact hourly, daily, monthly, or yearly units."
        )

    return freq