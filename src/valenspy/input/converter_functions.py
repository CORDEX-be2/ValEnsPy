"""A collection of functions unique to an input dataset required for conversion to valenspy compliant xarray datasets."""

from valenspy._utilities import load_yml, _fix_lat_lon
import xarray as xr

CORDEX_VARIABLES = load_yml("CORDEX_variables")


def EOBS_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw EOBS data to a ValensPy compliant xarray Dataset.

    Rename latitude and longitude coordinates to lat and lon, respectively.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of EOBS observations to convert

    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """
    ds = _fix_lat_lon(ds)

    return ds


def ERA5_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw ERA5 data to a ValensPy compliant xarray Dataset.

    Rename latitude and longitude coordinates to lat and lon, respectively. Rename valid_time to time for certain variables.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ERA5 observations to convert

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
    Convert a xarray with raw CCLM data to a ValensPy compliant xarray Dataset.

    Flatten the pressure dimension by renaming variables with the pressure level in the name and removing the pressure dimension.
    Drop the last time step of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CCLM simulation to convert

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


def ALARO_K_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw ALARO_K data to a ValensPy compliant xarray Dataset.

    Does nothing, WIP

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ALARO_K simulation to convert

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

def RADCLIM_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw RADCLIM data to a ValensPy compliant xarray Dataset.

    Rename nlon and nlat to lon and lat, respectively. Set the coordinates lat_bounds and lon_bounds as coordinates.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CCLM simulation to convert

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

def MAR_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw MAR data to a ValensPy compliant xarray Dataset.

    Rename TIME to time and remove the ZTQLEV and ZUVLEV dimensions by selecting the first value of each dimension.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of MAR simulation to convert
    
    Returns
    -------
    Dataset
        xarray dataset ready for unit conversion

    """
    ds = ds.rename({'TIME':'time'})
    ds = ds.isel(ZTQLEV=0,ZUVLEV=0)
    
    return ds