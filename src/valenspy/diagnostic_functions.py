import xarray as xr
from functools import wraps

###################################
# Model2Self diagnostic functions #
###################################


def diurnal_cycle(ds: xr.Dataset):
    """Calculate the diurnal cycle of the data. If lat and lon are present, the data is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the diurnal cycle of.

    Returns
    -------
    xr.Dataset
        The diurnal cycle of the data.
    """
    ds = _average_over_dims(ds, ["lat", "lon"])

    return ds.groupby("time.hour").mean("time")


def time_series_spatial_mean(ds: xr.Dataset):
    """Calculate the time series of the data. If lat and lon are present, the data is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the time series of the spatial mean of.

    Returns
    -------
    xr.Dataset
        The time series of the spatial mean of the data.
    """
    return _average_over_dims(ds, ["lat", "lon"])


##################################
# Model2Ref diagnostic functions #
##################################


def spatial_bias(ds: xr.Dataset, ref: xr.Dataset, calc_relative=False):
    """Calculate the spatial bias of the data compared to the reference. The time dimensions are averaged over if present.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the spatial bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    calc_relative : bool, optional
        If True, return the relative bias, if False return the absolute bias, by default False

    Returns
    -------
    xr.Dataset
        The spatial bias of the data compared to the reference.
    """
    return bias(
        _average_over_dims(ds, "time"),
        _average_over_dims(ref, "time"),
        calc_relative=calc_relative,
    )


def temporal_bias(ds: xr.Dataset, ref: xr.Dataset, calc_relative=False):
    """Calculate the temporal bias of the data compared to the reference. If lat and lon are present, ds and ref is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the temporal bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    calc_relative : bool, optional
        If True, return the relative bias, if False return the absolute bias, by default False

    Returns
    -------
    xr.Dataset
        The temporal bias of the data compared to the reference.
    """
    return bias(
        _average_over_dims(ds, ["lat", "lon"]),
        _average_over_dims(ref, ["lat", "lon"]),
        calc_relative=calc_relative,
    )


def diurnal_cycle_bias(ds: xr.Dataset, ref: xr.Dataset, calc_relative=False):
    """Calculate the diurnal cycle bias of the data compared to the reference. If lat and lon are present,  ds and ref is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the diurnal cycle bias of.
    ref : xr.Dataset
        The reference data to compare the data to.
    calc_relative : bool, optional
        If True, return the calc_relative bias, by default False

    Returns
    -------
    xr.Dataset
        The diurnal cycle bias of the data compared to the reference.
    """
    ds = _average_over_dims(ds, ["lat", "lon"])
    ref = _average_over_dims(ref, ["lat", "lon"])

    return bias(
        ds.groupby("time.hour").mean("time"),
        ref.groupby("time.hour").mean("time"),
        calc_relative=calc_relative,
    )


##################################
####### Helper functions #########
##################################


def _average_over_dims(da: xr.DataArray, dims):
    """Calculate the average over the specified dimensions if they are present in the data. Otherwise, return the data as is.

    Parameters
    ----------
    ds : xr.DataArray
        The data to calculate the spatial average of.
    dims : list or str
        The dimension(s) to average over.

    Returns
    -------
    xr.DataArray
        The data with the specified dimensions averaged over.
    """
    if isinstance(dims, str):
        dims = [dims]
    if all(dim not in ds.dims for dim in dims):
        return data
    return ds.mean([dim for dim in dims if dim in ds.dims], keep_attrs=True)


##################################
########### Metrics ##############
##################################


def bias(da: xr.DataArray, ref: xr.DataArray, calc_relative=False):
    """Calculate the bias of the data compared to a reference.

    Parameters
    ----------
    da : xr.DataArray
        The data to calculate the bias of.
    ref : xr.DataArray
        The reference to compare the data to.
    calc_relative : bool, optional
        If True, calculate the relative bias, if False calculate the absolute bias, by default False

    Returns
    -------
    xr.DataArray
        The bias of the data compared to there reference.
    """
    if calc_relative:
        return (da - ref) / ref
    else:
        return da - ref

######################################
############## Wrappers ##############
######################################

def requires_variables(variables):
    """
    A decorator that checks if the required variables are present in the dataset before applying the diagnostic.
    The required variables are specified as a list of strings. Only if all the required variables are present the diagnostic is applied.
    Note that this is a minimum requirement, the ds may contain more variables than the required ones.
    

    Parameters
    ----------
    variables : str or list of str
        The variable(s) required to apply the diagnostic.

    Example
    -----
    @requires_variables(["tas", "pr"])
    def my_diagnostic(ds: xr.Dataset):
        return ds.tas + ds.pr
    """
    def decorator(diagnostic_function):
        @wraps(diagnostic_function)
        def wrapper(ds, *args, **kwargs):
            if isinstance(variables, str):
                variables = [variables]
            if not all(var in ds for var in variables):
                raise ValueError(f"Variables {variables} are required to apply the diagnostic.")
            return diagnostic_function(ds, *args, **kwargs)
        return wrapper
    return decorator
