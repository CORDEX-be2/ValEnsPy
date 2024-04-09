import xarray as xr
##################################
# Model2Ref diagnostic functions #
##################################

def spatial_bias(data: xr.Dataset, ref: xr.Dataset, relative=False):
    """Calculate the spatial bias of the data compared to the reference.
    
    Parameters
    ----------
    data : DataTree
        The data to calculate the spatial bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    
    Returns
    -------
    Dataset
        The spatial bias of the data compared to the reference.
    """
    return bias(data.mean("time").tas, ref.mean("time").tas, relative)


##################################
########### Metrics ##############
##################################

def bias(data: xr.DataArray, ref: xr.DataArray, relative=False):
    """Calculate the bias of the data compared to a reference.
    
    Parameters
    ----------
    data : DataTree
        The data to calculate the bias of.
    ref : xr.Dataset
        The reference to compare the data to.
    
    Returns
    -------
    Dataset
        The bias of the data compared to there reference.
    """
    if relative:
        return (data - ref) / ref
    else:
        return data - ref

