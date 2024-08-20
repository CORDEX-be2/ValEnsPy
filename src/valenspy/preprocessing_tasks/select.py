# Collection of preprocessing functions to perfrom a selection

import xarray as xr
import numpy as np
from valenspy._regions import region_bounds

# make sure attributes are passed through
xr.set_options(keep_attrs=True)


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


def select_point(ds: xr.Dataset, lon_lat_point: tuple, rotated_pole: bool = False):
    """
    Select a point from the dataset based on the provided geographic coordinates.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset from which to select the point.
    lon_lat_point : tuple
        Geographic coordinates as a (longitude, latitude) pair in degrees.
    rotated_pole : bool, optional
        If True, the dataset is in a rotated pole projection and the input coordinates
        will be converted to rotated pole coordinates before selection. Default is False.

    Returns
    -------
    xarray.Dataset
        The dataset subset at the nearest point to the specified coordinates.
    """
    if rotated_pole:
        # Convert geographic coordinates to rotated pole coordinates
        lon_lat_point_rot = convert_geo_to_rot(lon_lat_point, ds)
        # Select the nearest point in the rotated pole coordinates
        ds_point = ds.sel(
            rlon=lon_lat_point_rot[0], rlat=lon_lat_point_rot[1], method="nearest"
        )
    else:
        # Select the nearest point based on geographic coordinates
        ds_point = ds.sel(lon=lon_lat_point[0], lat=lon_lat_point[1], method="nearest")

    return ds_point


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
    ds_region = select_region(ds, 'europe')
    """

    # get region bounds
    lat_bounds = region_bounds[region]["lat_bounds"]
    lon_bounds = region_bounds[region]["lon_bounds"]

    ds_sel = ds.sel(
        lon=slice(lon_bounds[0], lon_bounds[1]), lat=slice(lat_bounds[0], lat_bounds[1])
    )
    return ds_sel
