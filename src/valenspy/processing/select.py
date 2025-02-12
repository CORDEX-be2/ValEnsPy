# Collection of preprocessing functions to perfrom a selection

import xarray as xr
import numpy as np
import regionmask
import geopandas as gpd
from pathlib import Path
from valenspy._utilities._regions import region_bounds
import pyproj
import cartopy.crs as ccrs


# make sure attributes are passed through
xr.set_options(keep_attrs=True)


def select_region(ds: xr.Dataset, region: str):
    """
    Selects a specific geographical region from an xarray Dataset based on given region bounds.

    Parameters:
    ds (xr.Dataset): The input xarray Dataset from which to select the region.
    region (str): The name of the region to select. This should correspond to a key in the
                  `region_bounds` dictionary, which contains latitude and longitude bounds
                  for various regions.

    Returns:
    xr.Dataset: A new xarray Dataset containing only the data within the specified region.

    Example:
    ds_region = sel_region(ds, 'europe')
    """

    # get region bounds
    lat_bounds = region_bounds[region]["lat_bounds"]
    lon_bounds = region_bounds[region]["lon_bounds"]

    ds_sel = ds.sel(
        lon=slice(lon_bounds[0], lon_bounds[1]), lat=slice(lat_bounds[0], lat_bounds[1])
    )
    return ds_sel

def convert_geo_to_rot(coord: tuple, ds: xr.Dataset):
    """
    Converts a geographic (longitude, latitude) point to a rotated pole (rlon, rlat) point.

    Parameters
    ----------
    coord : list or tuple
        Geographic coordinates as a (longitude, latitude) pair in degrees.
    ds : xarray.Dataset
        The input dataset containing the rotated pole grid information (e.g., from a COSMO-CLM file).

    Returns
    -------
    list
        The corresponding rotated pole coordinates as [rlon, rlat] in degrees.
    """
    # Extract the rotated pole latitude and longitude from the dataset
    rp_lat = float(ds.rotated_pole.grid_north_pole_latitude)
    rp_lon = float(ds.rotated_pole.grid_north_pole_longitude)

    # Convert geographic and rotated pole coordinates from degrees to radians
    co = np.deg2rad(coord)
    rp_lat = np.deg2rad(rp_lat)
    rp_lon = np.deg2rad(rp_lon)

    # Calculate the rotated pole latitude (rlat) using spherical trigonometry
    p_rlat = np.arcsin(
        np.sin(co[1]) * np.sin(rp_lat)
        + np.cos(co[1]) * np.cos(rp_lat) * np.cos(co[0] - rp_lon)
    )

    # Calculate the rotated pole longitude (rlon) using spherical trigonometry
    p_rlon = np.arctan(
        (np.cos(co[1]) * np.sin(co[0] - rp_lon))
        / (
            np.cos(co[1]) * np.sin(rp_lat) * np.cos(co[0] - rp_lon)
            - np.sin(co[1]) * np.cos(rp_lat)
        )
    )

    # Convert the rotated pole coordinates from radians back to degrees
    p_rlat = np.rad2deg(p_rlat)
    p_rlon = np.rad2deg(p_rlon)

    # Return the rotated pole coordinates as a list
    return [p_rlon, p_rlat]

def convert_geo_to_LCC(coord: tuple, ds: xr.Dataset):
    """
    Convert the geographic coordinates to Lambert Conformal Coordinates.

    Parameters
    ----------
    coord : tuple
        Geographic coordinates as a (longitude, latitude) pair in degrees.
    ds : xarray.Dataset
        The input dataset containing the Lambert Conformal Conic grid information.
    
    Returns
    -------
    tuple
        The corresponding Lambert Conformal Coordinates as (x, y) in meters.
    """

    crs = ccrs.LambertConformal(
        central_longitude=ds.crs.longitude_of_central_meridian,
        central_latitude=ds.crs.latitude_of_projection_origin,
        false_easting=ds.crs.false_easting*1000,
        false_northing=ds.crs.false_northing*1000,
        standard_parallels=[ds.crs.standard_parallel]
    )

    transformer = pyproj.Transformer.from_crs(ccrs.PlateCarree(), crs)
    x,y = transformer.transform(*coord)
    return x/1000, y/1000

#TODO: fix this function to work using the ds.crs attribute so that not each crs has to be handled separately
#This adds the responsibility to the user to have a wel defined crs attribute (maybe some functionality to check this or help add this in input converter would be nice!)
def select_point(ds: xr.Dataset, lon_point: float, lat_point: float, projection: str = None):
    """
    Select a point from the dataset based on the provided geographic coordinates.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset from which to select the point.
    lon_point : float
        Geographic longitude coordinate 
    lat_point : float
        Geographic latitude coordinate   
    projection : str, optional
        The projection of the dataset. The point will be projected to the dataset's coordinate system before selection. 
        Currently supported projections are 'rotated_pole' and 'lcc' (Lambert Conformal Coordinates).

    Returns
    -------
    xarray.Dataset
        The dataset subset at the nearest point to the specified coordinates.
    """
    if projection == "rotated_pole":
        # Convert geographic coordinates to rotated pole coordinates
        lon_lat_point_rot = convert_geo_to_rot((lon_point,lat_point), ds)
        # Select the nearest point in the rotated pole coordinates
        ds_point = ds.sel(
            rlon=lon_lat_point_rot[0], rlat=lon_lat_point_rot[1], method="nearest"
        )
    elif projection == "lcc":
        point = convert_geo_to_LCC((lon_point,lat_point), ds)
        ds_point = ds.sel(x=point[0], y=point[1], method="nearest")
    else:
        # Select the nearest point based on geographic coordinates
        ds_point = ds.sel(lon=lon_point, lat=lat_point, method="nearest")

    return ds_point

