# collection of lat lon bounds of pre-defined regions - for plotting purposes in Valenspy

import xarray as xr 


# define region bounds bounds 
region_bounds ={
    'europe':
            {
            'lat_bounds' : [35, 70], 
            'lon_bounds' : [-15, 40]
            }, 
    'belgium': 
            {
            'lat_bounds' : [49, 52], 
            'lon_bounds' : [2, 7]
            }
    }


## is this the right location for this function?? Or should it go in preprocessing/or utilities? 
def sel_region(ds: xr.Dataset, region: str):
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
    lat_bounds = region_bounds[region]['lat_bounds']
    lon_bounds = region_bounds[region]['lon_bounds']

    ds_sel = ds.sel(lon=slice(lon_bounds[0],lon_bounds[1]),lat=slice(lat_bounds[0], lat_bounds[1]))
    return ds_sel