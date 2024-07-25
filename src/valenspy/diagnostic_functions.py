import xarray as xr

###################################
# Model2Self diagnostic functions #
###################################


def diurnal_cycle(data: xr.Dataset):
    """Calculate the diurnal cycle of the data. If lat and lon are present, the diurnal cycle is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    data : DataTree
        The data to calculate the diurnal cycle of.

    Returns
    -------
    Dataset
        The diurnal cycle of the data.
    """
    data = _average_over_dims(data, ["lat", "lon"])

    return data.groupby("time.hour").mean("time")


def time_series_spatial_mean(data: xr.Dataset):
    """Calculate the time series of the spatial mean of the data.

    Parameters
    ----------
    data : DataTree
        The data to calculate the time series of the spatial mean of.

    Returns
    -------
    Dataset
        The time series of the spatial mean of the data.
    """
    return _average_over_dims(data, ["lat", "lon"])


##################################
# Model2Ref diagnostic functions #
##################################


def spatial_bias(data: xr.Dataset, ref: xr.Dataset, compute_relative_bias=False):
    """Calculate the spatial bias of the data compared to the reference. Time dimensions are averaged over.

    Parameters
    ----------
    data : DataTree
        The data to calculate the spatial bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    compute_relative_bias : bool, optional
        If True, return the relative bias, if False return the absolute bias, by default False

    Returns
    -------
    Dataset
        The spatial bias of the data compared to the reference.
    """
    return bias(
        _average_over_dims(data, "time"),
        _average_over_dims(ref, "time"),
        compute_relative_bias=compute_relative_bias,
    )


def temporal_bias(data: xr.Dataset, ref: xr.Dataset, compute_relative_bias=False):
    """Calculate the temporal bias of the data compared to the reference. Spatial dimensions are averaged over.

    Parameters
    ----------
    data : DataTree
        The data to calculate the temporal bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    compute_relative_bias : bool, optional
        If True, return the relative bias, if False return the absolute bias, by default False

    Returns
    -------
    Dataset
        The temporal bias of the data compared to the reference.
    """
    return bias(
        _average_over_dims(data, ["lat", "lon"]),
        _average_over_dims(ref, ["lat", "lon"]),
        compute_relative_bias=compute_relative_bias,
    )


def diurnal_cycle_bias(data: xr.Dataset, ref: xr.Dataset, compute_relative_bias=False):
    """Calculate the diurnal cycle bias of the data compared to the reference. If lat and lon are present, the diurnal cycle is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    data : DataTree
        The data to calculate the diurnal cycle bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    compute_relative_bias : bool, optional
        If True, return the compute_relative_bias bias, by default False

    Returns
    -------
    Dataset
        The diurnal cycle bias of the data compared to the reference.
    """
    data = _average_over_dims(data, ["lat", "lon"])
    ref = _average_over_dims(ref, ["lat", "lon"])

    return data.groupby("time.hour").mean("time") - ref.groupby("time.hour").mean(
        "time"
    )


##################################
####### Helper functions #########
##################################


def _average_over_dims(data: xr.Dataset, dims):
    """Calculate the average over the specified dimensions if they are present in the data. Otherwise, return the data as is.

    Parameters
    ----------
    data : DataTree
        The data to calculate the spatial average of.
    dims : list or str
        The dimension(s) to average over.

    Returns
    -------
    Dataset
        The data with the specified dimensions averaged over.
    """
    if isinstance(dims, str):
        dims = [dims]
    if all(dim not in data.dims for dim in dims):
        return data
    return data.mean([dim for dim in dims if dim in data.dims], keep_attrs=True)


##################################
########### Metrics ##############
##################################


def bias(data: xr.DataArray, ref: xr.DataArray, compute_relative_bias=False):
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
    if compute_relative_bias:
        return (data - ref) / ref
    else:
        return data - ref
