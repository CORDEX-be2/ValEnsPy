import xarray as xr
import numpy as np


########### To Kelvin ############
def convert_Celcius_to_Kelvin(da: xr.DataArray):
    """
    Convert values in xarray DataArray from °C to K
    """
    da = da + 273.15  # Celcius to Kelvin
    da.attrs["units"] = "K"  # Naming as in CORDEX_variables.yml

    return da


########### To Celsius ############
def convert_Kelvin_to_Celcius(da: xr.DataArray):
    """
    Convert values in xarray DataArray from K to °C
    """
    da = da - 273.15  # Kelvin to Celcius
    da.attrs["units"] = "°C"

    return da


########### To Pa ############
def convert_hPa_to_Pa(da: xr.DataArray):
    """
    Convert values in xarray DataArray from hPa to Pa
    """
    da = da * 100
    da.attrs["units"] = "Pa"

    return da


########### To hPa ############
def convert_Pa_to_hPa(da: xr.DataArray):
    """
    Convert values in xarray DataArray from Pa to hPa
    """
    da = da / 100
    da.attrs["units"] = "hPa"

    return da


########### To kg m-2 s-1 ############
def convert_mm_to_kg_m2s(da: xr.DataArray):
    """
    Convert daily (!) values in xarray DataArray from mm to kg m^-2 s^-1
    """
    # first, get timestep (frequency) by calculating the difference between the first consecutive time values in seconds
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    da = da / timestep_nseconds  # mm to kg m^-2 s^-1
    da.attrs["units"] = "kg m-2 s-1"

    return da


def convert_mm_hr_to_kg_m2s(da: xr.DataArray):
    """
    Convert daily (!) values in xarray DataArray from mm to kg m^-2 s^-1
    """
    da = da / 3600  # mm to kg m^-2 s^-1
    da.attrs["units"] = "kg m-2 s-1"

    return da


def convert_kg_m2_to_kg_m2s(da: xr.DataArray):
    """
    Convert daily (!) values in xarray DataArray from kg m^-2 to kg m^-2 s^-1
    """
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    da = da / timestep_nseconds  # kg m^-2 to kg m^-2 s^-1
    da.attrs["units"] = "kg m-2 s-1"

    return da


def convert_m_to_kg_m2s(da: xr.DataArray):
    """
    Convert values in xarray DataArray from mm hr^-1 to kg m^-2 s^-1
    """
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")
    da = da * 1000 / timestep_nseconds  # mm to kg m^-2 s^-1
    da.attrs["units"] = "kg m-2 s-1"

    return da


def convert_m_hr_to_kg_m2s(da: xr.DataArray):
    """
    Convert values in xarray DataArray from m hr^-1 to kg m^-2 s^-1s
    """
    # do conversion
    da = da * 1000 / 3600  # m hr^-1 to kg m^-2 s^-1
    da.attrs["units"] = "kg m-2 s-1"

    return da


########### To m h-1 ############
def convert_kg_m2s_to_mh(da: xr.DataArray):
    """
    Convert values in xarray DataArray from kg m^-2 s^-1 to mm hr^-1
    """
    da = da * 3600 / 1000  # kg m^-2 s^-1 to mm hr^-1
    da.attrs["units"] = "mm hr-1"

    return da


########### To W m-2 ############
def convert_J_m2_to_W_m2(da: xr.DataArray):
    """
    Convert values in xarray DataArray from J m^2 to W m^2
    """
    # first, get timestep (frequency) by calculating the difference between the first consecutive time values in seconds
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    da = da / timestep_nseconds  # J m^2 to W m^2
    da.attrs["units"] = "W m-2"

    return da


def convert_kWh_m2_day_to_W_m2(da: xr.DataArray):
    """
    Convert values in xarray DataArray from kWh/m2/day to W m^2
    """
    da = da * (1000) / 24
    da.attrs["units"] = "W m-2"

    return da


########### To % ############
def convert_fraction_to_percent(da: xr.DataArray):
    """
    Convert values in xarray DataArray from unitless to %
    """
    da = da * 100
    da.attrs["units"] = "%"

    return da


########### Helper functions ############
def _determine_time_interval(da: xr.DataArray):
    """
    Find the time interval (freq) of the input data array based on it's time axis, by calculating the difference between the first two time instances.
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
