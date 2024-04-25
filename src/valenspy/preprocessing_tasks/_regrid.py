from valenspy.preprocessing_tasks.task import PreprocessingTask
import xarray as xr
import numpy as np
from typing import Union


class Regrid(PreprocessingTask):
    REGRID_METHODS_XARRAY = {"linear", "nearest"} #Native xarray regridding methods (for n-dimensional data)
    """A regrid preprocessing task."""

    def __init__(self, target_grid: Union[xr.Dataset, str], name="", description=None):
        """Initialize the Regrid task.
        
        Parameters
        ----------
        target_grid : xr.Dataset or str
            The target grid to regrid to. This can be a string representation of a global grid resolution (e.g. "1x1") 
            or an xr.Dataset with lat and lon coordinates.
        """
        super().__init__("regrid_" + name, description)

        if isinstance(target_grid, str):
            self.target_grid = create_global_grid(*parse_global_grid(target_grid))
        else:
            self.target_grid = target_grid

    def apply(self, data: xr.Dataset, method="linear") -> xr.Dataset:
        """Apply the regridding task to an xr.Dataset.

        Parameters
        ----------
        data : xr.Dataset
            The data to regrid from.

        Returns
        -------
        xr.Dataset
            The regridded dataset.
        """

        if method not in self.REGRID_METHODS_XARRAY:
            raise NotImplementedError(
                f"{method} regridding is not yet implemented. Please choose from {self.REGRID_METHODS_XARRAY}."
            )

        return data.interp(
            lon=self.target_grid.lon, lat=self.target_grid.lat, method=method
        )


def parse_global_grid(mxn="1x1"):
    """Parse a string representation of a global grid resolution.

    Parameters
    ----------
    mxn : str, optional
        The string representation of the grid resolution, by default "1x1".

    Returns
    -------
    tuple
        The grid resolution as a tuple of integers.
    """
    m, n = mxn.split("x")
    return int(m), int(n)

# Stock xarray - global grid extents (degrees).
_LAT_MIN = -90.0
_LAT_MAX = 90.0
_LAT_RANGE = _LAT_MAX - _LAT_MIN
_LON_MIN = 0.0
_LON_MAX = 360.0
_LON_RANGE = _LON_MAX - _LON_MIN

def create_global_grid(m, n, lat_offset=True, lon_offset=True):
    """Create a global lat-lon grid on a standard mxn grid.

    The longitude range is from 0 to 360 degrees, and the latitude range is from -90 to 90 degrees. 

    Parameters
    ----------
    m : int
        The latitudinal resolution.
    n : int
        The longitudinal resolution.
    lat_offset : bool, optional
        Whether to offset the latitude by m degrees, by default True.
    lon_offset : bool, optional
        Whether to offset the longitude by n degrees, by default True.

    Returns
    -------
    xr.Dataset
        The global grid.
    """
    lat_offset_val = m / 2 if lat_offset else 0
    lon_offset_val = n / 2 if lon_offset else 0

    lat = np.linspace(_LAT_MIN + lat_offset_val, _LAT_MAX - lat_offset_val, num=int(_LAT_RANGE//m))
    lon = np.linspace(_LON_MIN + lon_offset_val, _LON_MAX - lon_offset_val, num=int(_LON_RANGE//n))

    grid = xr.Dataset(
        coords={
            "lat": lat,
            "lon": lon,
        }
    )
    grid["lat"].attrs = {"units": "degrees_north", "long_name": "latitude"}
    grid["lon"].attrs = {"units": "degrees_east", "long_name": "longitude"}

    #TODO: Add lat_bnds and lon_bnds

    return grid


