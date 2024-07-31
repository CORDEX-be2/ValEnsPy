import xarray as xr

# make ure attributes are passed through
xr.set_options(keep_attrs=True)  

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

from functools import wraps

###################################
# Model2Self diagnostic functions #
###################################


def diurnal_cycle(ds: xr.Dataset):
    """Calculate the diurnal cycle of the data. If lat and lon are present, the data is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset or xr.DataArray
        The data to calculate the diurnal cycle of.

    Returns
    -------
    xr.Dataset or xr.DataArray
        The diurnal cycle of the data.
    """
    ds = _average_over_dims(ds, ["lat", "lon"])

    return ds.groupby("time.hour").mean("time")


def time_series_spatial_mean(ds: xr.Dataset):
    """Calculate the time series of the data. If lat and lon are present, the data is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset or xr.DataArray
        The data to calculate the time series of the spatial mean of.

    Returns
    -------
    xr.Dataset or xr.DataArray
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
    ds : xr.Dataset or xr.DataArray
        The data to calculate the spatial bias of.
    ref : xr.Dataset or xr.DataArray
        The reference data to compare the data to.
    calc_relative : bool, optional
        If True, return the relative bias, if False return the absolute bias, by default False

    Returns
    -------
    xr.Dataset or xr.DataArray
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
    ds : xr.Dataset or xr.DataArray
        The data to calculate the temporal bias of.
    ref : xr.Dataset or xr.DataArray
        The reference data to compare the data to.
    calc_relative : bool, optional
        If True, return the relative bias, if False return the absolute bias, by default False

    Returns
    -------
    xr.Dataset or xr.DataArrays
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
    ds : xr.Dataset or xr.DataArray
        The data to calculate the diurnal cycle bias of.
    ref : xr.Dataset or xr.DataArray
        The reference data to compare the data to.
    calc_relative : bool, optional
        If True, return the calc_relative bias, by default False

    Returns
    -------
    xr.Dataset or xr.DataArray
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
    ds : xr.DataArray or xr.Dataset
        The data to calculate the spatial average of.
    dims : list or str
        The dimension(s) to average over.

    Returns
    -------
    xr.DataArray or xr.Dataset
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
    da : xr.DataArray or xr.Dataset
        The data to calculate the bias of.
    ref : xr.DataArray or xr.Dataset
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
    A decorator that checks if the required variables are present in the dataset (and reference dataset if applicable) before applying the diagnostic.
    The required variables are specified as a list of strings. Only if all the required variables are present the diagnostic is applied.
    Note that this is a minimum requirement, the ds may contain other variables than the required ones.

    Parameters
    ----------
    variables : str or list of str
        The variable(s) required to apply the diagnostic.

    Examples
    -----
    #The diagnostic function requires the variables 'tas' and 'pr' to be present in the dataset.
    @requires_variables(["tas", "pr"])
    def my_diagnostic(ds: xr.Dataset):
        return ds.tas + ds.pr

    #This also checks if the variables are present in both the data and the reference.
    #An error is raised if the required variables are not present in the data or the reference.
    @requires_variables(["tas", "pr"])
    def my_diagnostic(ds: xr.Dataset, ref: xr.Dataset):
        return ds.tas + ref.pr
    """

    def decorator(diagnostic_function):
        @wraps(diagnostic_function)
        def wrapper(ds, *args, **kwargs):
            required_vars = [variables] if isinstance(variables, str) else variables
            # Do the check for the ds
            if not all(var in ds.variables for var in required_vars):
                raise ValueError(
                    f"Variables {required_vars} are required to apply the diagnostic."
                )
            # Do the check for the reference if it is present, the reference is the second argument after the ds argument and should be a xr.Dataset.
            if len(args) > 0 and isinstance(args[0], xr.Dataset):
                ref = args[0]
                if not all(var in ref.variables for var in required_vars):
                    raise ValueError(
                        f"Variables {required_vars} are required to apply the diagnostic."
                    )
            return diagnostic_function(*args, **kwargs)

        return wrapper

    return decorator
