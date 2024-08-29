# collection of lat lon bounds of pre-defined regions - for plotting purposes in Valenspy

import xarray as xr


# define region bounds bounds
region_bounds = {
    "europe": {"lat_bounds": [35, 70], "lon_bounds": [-15, 40]},
    "belgium": {"lat_bounds": [49, 52], "lon_bounds": [2, 7]},
}
