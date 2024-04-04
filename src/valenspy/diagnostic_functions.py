import xarray as xr
##################################
# Model2Obs diagnostic functions #
##################################

def spatial_bias(data: xr.Dataset, obs: xr.Dataset, relative=False):
    """Calculate the spatial bias of the data compared to the observations.
    
    Parameters
    ----------
    data : DataTree
        The data to calculate the spatial bias of.
    obs : xr.Dataset
        The observations to compare the data to.
    
    Returns
    -------
    Dataset
        The spatial bias of the data compared to the observations.
    """
    return bias(data.mean("time").tas, obs.mean("time").tas, relative)




##################################
########### Metrics ##############
##################################

def bias(data: xr.DataArray, obs: xr.DataArray, relative=False):
    """Calculate the bias of the data compared to the observations.
    
    Parameters
    ----------
    data : DataTree
        The data to calculate the bias of.
    obs : xr.Dataset
        The observations to compare the data to.
    
    Returns
    -------
    Dataset
        The bias of the data compared to the observations.
    """
    if relative:
        return (data - obs) / obs
    else:
        return data - obs

