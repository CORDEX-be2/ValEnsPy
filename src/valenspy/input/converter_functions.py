from pathlib import Path
from valenspy._utilities import load_yml, is_cf_compliant, cf_status
from valenspy.input.unit_converter import convert_all_units_to_CF, _determine_time_interval
import xarray as xr
import pandas as pd
import numpy as np

CORDEX_VARIABLES = load_yml("CORDEX_variables")

def _set_global_attributes(ds: xr.Dataset, metadata_info):
    for key, value in metadata_info.items():
        ds.attrs[key] = value

    return ds

def _fix_lat_lon(ds: xr.Dataset):
    # rename dimensions if not yet renamed
    if "lon" not in ds.coords:
        ds = ds.rename({"longitude": "lon"})
    if "lat" not in ds.coords:
        ds = ds.rename({"latitude": "lat"})

    # make sure lat and lon are sorted ascending
    ds = ds.sortby("lat").sortby("lon")
    return ds


def EOBS_to_CF(ds: xr.Dataset) -> xr.Dataset:
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
    ds = _fix_lat_lon(ds)

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

    # bugfix ERA5 (found in clh): replace valid_time by time
    if "time" not in ds:
        ds = ds.rename({"valid_time": "time"})

    ds = _fix_lat_lon(ds)

    return ds


def CCLM_to_CF(ds: xr.Dataset) -> xr.Dataset:
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

    #For each variable in the dataset which has a pressure dimension, create a new variable with the pressure level in the name and remove the pressure dimension
    if "pressure" in ds.dims:
        for var in ds.data_vars:
            if "pressure" in ds[var].dims:
                for pressure in ds[var].pressure.values:
                    new_var = var + str(int(pressure / 100)) + "p"
                    ds[new_var] = ds[var].sel(pressure=pressure)
                ds = ds.drop_vars(var)
        ds = ds.drop_dims("pressure")

    # Seems to be failing because the last time step of the dataset is at 11:00:00, while all others are at 11:30:00.
    # One option: 
    ds = ds.isel(time=slice(0,-1))
    # Other options:
    # ds = ds.assign_coords(time=(ds.time.dt.floor('H'))) Note ciel does not work as the last time step is at 11:00:00 and the ceil function will round it to 12:00:00
    # Original option (only for hourly data)
    # new_time = ds_cclm.tas.time.astype('datetime64[D]') + np.timedelta64(12, 'h')
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
    # # Special conversion for precipitation
    # if "rain_convective" in ds.data_vars and "rain_stratiform" in ds.data_vars:
    #     ds["pr"] = ds["rain_convective"] + ds["rain_stratiform"]
    #     ds["pr"].attrs["units"] = "mm"
    #     ds["pr"] = convert_mm_to_kg_m2s(ds["rain_convective"] + ds["rain_stratiform"])
    #     ds["pr"].attrs["standard_name"] = "precipitation_flux"
    #     ds["pr"].attrs["long_name"] = "Precipitation"
    #     ds["pr"].attrs["dataset"] = model_name
    #     ds["pr"].attrs["original_name"] = "rain_convective + rain_stratiform"
    #     for key, value in metadata_info.items():
    #         ds["pr"].attrs[key] = value

    #     # Assuming monthly decumilation! This is not always the case!
    #     def decumilate(ds):
    #         ds_decum = ds.diff("time")
    #         # Add the first value of the month of original dataset to the decumilated dataset
    #         ds_decum = xr.concat([ds.isel(time=0), ds_decum], dim="time")
    #         return ds_decum

    #     ds.coords["year_month"] = ds["time.year"] * 100 + ds["time.month"]
    #     ds["pr"] = ds["pr"].groupby("year_month").apply(decumilate)

    return ds

def RADCLIM_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the RADCLIM xarray netCDF to a CF compliant xarray Dataset

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

    ds = ds.set_coords(("lat_bounds", "lon_bounds"))

    if "nlon" in ds.dims:
        ds = ds.rename({"nlon": "lon"})
    if "nlat" in ds.dims:
        ds = ds.rename({"nlat": "lat"})

    return ds
