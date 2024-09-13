import xarray as xr
import numpy as np
from scipy.stats import spearmanr


# make sure attributes are passed through
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


##################################
# Model2Ref diagnostic functions #
##################################


def spatial_bias(ds: xr.Dataset, ref: xr.Dataset, calc_relative=False):
    """Calculate the spatial bias of the data compared to the reference. The time dimensions are averaged over if present.

    Parameters
    ----------
    ds : xr.Dataset
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


def _average_over_dims(ds: xr.Dataset, dims):
    """Calculate the average over the specified dimensions if they are present in the data. Otherwise, return the data as is.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the spatial average of.
    dims : list or str
        The dimension(s) to average over.

    Returns
    -------
    xr.Dataset
        The data with the specified dimensions averaged over.
    """
    if isinstance(dims, str):
        dims = [dims]
    if all(dim not in ds.dims for dim in dims):
        return ds
    return ds.mean([dim for dim in dims if dim in ds.dims], keep_attrs=True)


################################################
########### Metrics & skill scores #############
################################################


def bias(da: xr.Dataset, ref: xr.Dataset, calc_relative=False):
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
    xr.Datasets
        The bias of the data compared to there reference.
    """
    if calc_relative:
        return (da - ref) / ref
    else:
        return da - ref

def mean_bias(da_mod: xr.Dataset, da_ref: xr.Dataset):
    """Calculate the bias of the means of modeled and reference data.

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
    xr.Datasets
        The bias of the data compared to there reference.
    """
    return (da_mod - da_ref).mean().values

def mean_absolute_error(da_mod: xr.DataArray, da_ref: xr.DataArray, percentile: float = None) -> float:
    """
    Calculate the Mean Absolute Error (MAE) between model forecasts and reference data.
    Optionally, calculate the MAE based on a specified percentile.

    Parameters
    ----------
    da_mod : xr.DataArray
        The model forecast data to compare.
    da_ref : xr.DataArray
        The reference data to compare against.
    percentile : float, optional
        The percentile (0 to 1) to calculate the MAE for, using the quantile values of the data arrays.
        If None, calculates the MAE for the entire data without considering percentiles.

    Returns
    -------
    float
        The Mean Absolute Error (MAE) between the model and reference data, or at the specified percentile.
    """
    # Ensure the DataArrays have the same shape
    if da_mod.shape != da_ref.shape:
        raise ValueError("Model and reference data must have the same shape.")

    if percentile is None:
        # Calculate the MAE for the entire data
        mae = np.mean(np.abs(da_mod.values - da_ref.values))
    else:
        # Calculate the MAE for the specified percentile
        mod_percentile = da_mod.quantile(percentile)
        ref_percentile = da_ref.quantile(percentile)
        mae = np.mean(np.abs(mod_percentile.values - ref_percentile.values))
    
    return mae

def root_mean_square_error(da_mod: xr.DataArray, da_ref: xr.DataArray) -> float:
    """
    Calculate the Root Mean Square Error (RMSE) between model data and reference data.

    Parameters
    ----------
    da_mod : xr.DataArray
        The model data to compare (should match the shape of da_ref).
    da_ref : xr.DataArray
        The reference data to compare against (should match the shape of da_mod).

    Returns
    -------
    float
        The Root Mean Square Error (RMSE) between the model and reference data.
    """
    # Ensure the DataArrays have the same shape
    if da_mod.shape != da_ref.shape:
        raise ValueError("Model and reference data must have the same shape.")

    # Calculate the squared differences
    squared_diff = (da_mod - da_ref) ** 2
    
    # Calculate the mean of squared differences
    mean_squared_diff = squared_diff.mean().values
    
    # Calculate and return the RMSE
    rmse = np.sqrt(mean_squared_diff)
    
    return rmse


def spearman_correlation(da_mod: xr.DataArray, da_ref: xr.DataArray) -> float:
    """
    Calculate Spearman's rank correlation coefficient between model data and reference data.

    Parameters
    ----------
    da_mod : xr.DataArray
        The model data to compare (2D array where rows are observations and columns are variables).
    da_ref : xr.DataArray
        The reference data to compare (2D array where rows are observations and columns are variables).

    Returns
    -------
    float
        Spearman's rank correlation coefficient between the flattened model and reference data.
    """
    # Flatten the DataArrays to 1D arrays for correlation calculation
    mod_data = da_mod.values.flatten()
    ref_data = da_ref.values.flatten()
    
    # Ensure that the flattened arrays have the same length
    if len(mod_data) != len(ref_data):
        raise ValueError("Model and reference data must have the same length after flattening.")

    # Calculate Spearman's rank correlation
    correlation, _ = spearmanr(mod_data, ref_data)
    
    return correlation


def optimal_bin_width(da_mod: xr.DataArray, da_ref: xr.DataArray) -> float:
    """
    Calculate the optimal bin width for both forecast (da_mod) and observed (da_ref) data.
    
    Parameters:
    da_mod (xr.DataArray): Forecasted temperatures (continuous).
    da_ref (xr.DataArray): Observed temperatures (continuous).
    
    Returns:
    float: Optimal bin width for both datasets.
    """
    
    # Combine both datasets
    combined_data = xr.concat([da_mod, da_ref], dim="time")

    # Freedman-Diaconis rule: Bin width = 2 * (IQR / n^(1/3))
    q25 = combined_data.quantile(0.25).item()
    q75 = combined_data.quantile(0.75).item()
    iqr = q75 - q25
    n = combined_data.size
    bin_width = 2 * (iqr / np.cbrt(n))

    std_dev = np.std(combined_data)
    n = len(combined_data)
    binwdth = 3.5 * (std_dev / np.cbrt(n))
    return bin_width
    
   # return bin_width


def perkins_skill_score(da: xr.DataArray, ref: xr.DataArray, binwidth: float = None):
    """
    Calculate the Perkins Skill Score (PSS).

    Parameters
    ----------
    da : xr.DataArray
        The model data to compare.
    ref : xr.DataArray
        The reference data to compare against.
    binwidth : float
     The width of each bin for the histogram. If not provided, it is calculated.

    Returns
    -------
    float
        The Perkins Skill Score (PSS).
    """
    if binwidth is None: 
       binwidth  = optimal_bin_width(da, ref)

    # Flatten the DataArrays to 1D for comparison
    mod_data = da.values.flatten()
    ref_data = ref.values.flatten()

    # Define the edges of the bins based on the data range
    lower_edge = min(np.min(mod_data), np.min(ref_data))
    upper_edge = max(np.max(mod_data), np.max(ref_data))

    # Calculate the histograms
    freq_m, _ = np.histogram(mod_data, bins=np.arange(lower_edge, upper_edge + binwidth, binwidth))
    freq_r, _ = np.histogram(ref_data, bins=np.arange(lower_edge, upper_edge + binwidth, binwidth))

    # Normalize the histograms by the total number of data points to compare probabilities
    freq_m = freq_m / np.sum(freq_m)
    freq_r = freq_r / np.sum(freq_r)
    # Calculate and return the PSS
    return np.sum(np.minimum(freq_m, freq_r)), freq_m, freq_r, binwidth


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
