from pathlib import Path
from valenspy.cf_checks import is_cf_compliant, cf_status
from valenspy._utilities import load_yml
from valenspy._unit_conversions import convert_all_units_to_CF, _determine_time_interval, _convert_mm_to_kg_m2s
import xarray as xr
import pandas as pd
import numpy as np

CORDEX_VARIABLES = load_yml("CORDEX_variables")

def _set_global_attributes(ds: xr.Dataset, metadata_info):
    for key, value in metadata_info.items():
        ds.attrs[key] = value
    
    return ds

def _rename_to_lon(ds: xr.Dataset, lon_name):
    return ds.rename_vars({ds.lon.name: lon_name})

def _rename_to_lat(ds: xr.Dataset, lat_name):
    return ds.rename_vars({ds.lat.name: lat_name})

def EOBS_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the EOBS xarray netCDF to a CF compliant xarray Dataset

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of EOBS observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """
    #Unique information regarding the EOBS dataset
    obsdata_name = "EOBS"
    raw_LOOKUP = load_yml(f"{obsdata_name}_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq":"daily", "spatial_resolution":"0.1deg", "region":"europe"}

    #Convert all units to CF, add metadata and set global attributes
    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)
    
    # rename dimensions if not yet renamed
    if "lon" not in ds.coords:
        ds = ds.rename({"longitude": "lon"})
    if "lat" not in ds.coords:
        ds = ds.rename({"latitude": "lat"})

    # make sure lat and lon are sorted ascending
    ds = ds.sortby('lat').sortby('lon')

    # Soft check for CF compliance
    cf_status(ds)

    return ds

def ERA5_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the ERA5 xarray dataset to a xarray Dataset in CF convention

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ERA5 observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant ERA5 observations for the specified variable.
    """
    #Unique information regarding the ERA5 dataset
    obsdata_name = "ERA5"
    raw_LOOKUP = load_yml(f"{obsdata_name}_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq": _determine_time_interval(ds)}

    #Convert all units to CF, add metadata and set global attributes
    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)

    # rename dimensions if not yet renamed
    if "lon" not in ds.coords:
        ds = ds.rename({"longitude": "lon"})
    if "lat" not in ds.coords:
        ds = ds.rename({"latitude": "lat"})

    # make sure lat and lon are sorted ascending
    ds = ds.sortby('lat').sortby('lon')

    cf_status(ds)

    return ds

def ERA5Land_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the ERA5-Land xarray dataset to a xarray Dataset in CF convention.
    Uses the same lookup table as ERA5.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ERA5-Land observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant ERA5-Land observations for the specified variable.
    """
    #Unique information regarding the ERA5 dataset
    obsdata_name = "ERA5-Land"
    raw_LOOKUP = load_yml(f"ERA5_lookup") #Note that the ERA5 lookup is used for both ERA5 and ERA5-Land

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq": _determine_time_interval(ds)}

    #Convert all units to CF, add metadata and set global attributes
    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)

    # rename dimensions if not yet renamed
    if "lon" not in ds.coords:
        ds = ds.rename({"longitude": "lon"})
    if "lat" not in ds.coords:
        ds = ds.rename({"latitude": "lat"})
    
    # make sure lat and lon are sorted ascending
    ds = ds.sortby('lat').sortby('lon')    
    
    cf_status(ds)

    return ds

def CLIMATE_GRID_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the CLIMATE_GRID xarray dataset to a xarray Dataset in CF convention

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CLIMATE_GRID observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant CLIMATE_GRID observations for the specified variable.
    """
    obsdata_name = "CLIMATE_GRID"
    raw_LOOKUP = load_yml(f"{obsdata_name}_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq":"daily", "spatial_resolution":"0.07° x 0.045° (~5km)", "region":"belgium"}

    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)
    
    cf_status(ds)

    return ds

def CCLM_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the CCLM xarray netCDF to a CF compliant xarray Dataset

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CCLM simulation to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator?
    model_name = "CCLM"
    src_path = Path("../src/valenspy")

    # open observational specific lookup dictionary
    raw_LOOKUP = load_yml(f"{model_name}_lookup")

    if metadata_info is None: # Set standard metadata if not provided
        metadata_info = {"experiment":""}

    metadata_info["dataset"] = model_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)

    # set attributes in whole dataset
    ds = _set_global_attributes(ds, metadata_info)

    # Soft check for CF compliance
    cf_status(ds)

    return ds

def ALARO_K_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the ALARO_K (converted with Kwintens R scripts) from xarray netCDF to a CF compliant xarray Dataset

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ALARO_K simulation to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """
    model_name = "ALARO-SFX_K"
    raw_LOOKUP = load_yml(model_name + "_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {}

    metadata_info["dataset"] = model_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)

    #Special conversion for precipitation
    if "rain_convective" in ds.data_vars and "rain_stratiform" in ds.data_vars:
        ds["pr"] = _convert_mm_to_kg_m2s(ds["rain_convective"] + ds["rain_stratiform"])
        ds["pr"].attrs["standard_name"] = "precipitation_flux"
        ds["pr"].attrs["long_name"] = "Precipitation"
        ds["pr"].attrs["dataset"] = model_name
        ds["pr"].attrs["original_name"] = "rain_convective + rain_stratiform"
        for key, value in metadata_info.items():
            ds["pr"].attrs[key] = value

    ds = _set_global_attributes(ds, metadata_info)
    
    cf_status(ds)

    return ds

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