import xarray as xr
import numpy as np
from scipy.stats import spearmanr
from datatree import DataTree
import pandas as pd
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
        mae = np.nanmean(np.abs(da_mod.values - da_ref.values))
    else:
        # Calculate the MAE for the specified percentile
        mod_percentile = da_mod.compute().quantile(percentile)
        ref_percentile = da_ref.compute().quantile(percentile)
        mae = np.nanmean(np.abs(mod_percentile.values - ref_percentile.values))
    
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
    correlation, _ = spearmanr(mod_data, ref_data, nan_policy='omit')
    
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
    combined_data = xr.concat([da_mod, da_ref], dim="time").compute()

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
    
def get_userdefined_binwidth(variable):
    """
    Get user defined, hard coded binwidths for Perkins Skill score calculation
    """
    # define bin width lookup table
    d_binwidth = { 
    'tas'    : 2,
    'tasmax' : 2,
    'tasmin' : 2,
    'ps'     : 500,
    'psl'    : 500,
    'clt'    : 10,
    'clh'    : 10,
    'clm'    : 10,
    'cll'    : 10 }

    if variable in d_binwidth: 
        return d_binwidth[variable]
    else:

        print(f"{variable} has no defined binwidths")



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
    lower_edge = min(np.nanmin(mod_data), np.nanmin(ref_data))
    upper_edge = max(np.nanmin(mod_data), np.nanmin(ref_data))

    # Calculate the histograms
    freq_m, _ = np.histogram(mod_data, bins=np.arange(lower_edge, upper_edge + binwidth, binwidth))
    freq_r, _ = np.histogram(ref_data, bins=np.arange(lower_edge, upper_edge + binwidth, binwidth))

    # Normalize the histograms by the total number of data points to compare probabilities
    freq_m = freq_m / np.sum(freq_m)
    freq_r = freq_r / np.sum(freq_r)
    # Calculate and return the PSS
    return np.sum(np.minimum(freq_m, freq_r)), freq_m, freq_r, binwidth

def calc_metrics(dt_mod: DataTree, da_obs: xr.DataArray, save_csv=False, csv_path=None):
    """
    Calculate statistical performance metrics for model data against observed data.

    This function computes various metrics between the model data stored in the DataTree 
    object `dt_mod` and the observed data `da_obs`. Metrics include Mean Bias, Mean Absolute 
    Error (MAE) at different percentiles, Root Mean Square Error (RMSE), Spearman Correlation, 
    and Perkins Skill Score (PSS). The metrics are collected for each member of the `dt_mod` 
    tree, and a merged pandas DataFrame containing these metrics is returned. Optionally, the 
    resulting DataFrame can be saved as a CSV file.

    Parameters:
    -----------
    dt_mod : DataTree
        A DataTree containing the model data for different members. The function loops 
        through each member to calculate the metrics.
    da_obs : xr.DataArray
        The observed data to compare against the model data.
    save_csv : bool, optional
        If True, saves the resulting DataFrame as a CSV file. Default is False.
    csv_path : str, optional
        The path where the CSV file will be saved. If None, a default path is generated 
        using the variable names.

    Returns:
    --------
    df_metric : pd.DataFrame
        A DataFrame containing the calculated metrics for each member in the `dt_mod`.

    Metrics:
    --------
    - Mean Bias
    - Mean Absolute Error
    - MAE at 90th Percentile
    - MAE at 99th Percentile
    - MAE at 10th Percentile
    - MAE at 1st Percentile
    - Root Mean Square Error
    - Spearman Correlation
    - Perkins Skill Score
    """
    
    for i, member in enumerate(dt_mod):
        ds_mod = dt_mod[member].ds
        variable = list(ds_mod.keys())[0]
        da_mod = ds_mod[variable]
        
        # Calculate metrics
        bias = mean_bias(da_mod, da_obs)
        mae = mean_absolute_error(da_mod, da_obs)
        mae_90pctl = mean_absolute_error(da_mod, da_obs, percentile=0.9)
        mae_99pctl = mean_absolute_error(da_mod, da_obs, percentile=0.99)
        mae_10pctl = mean_absolute_error(da_mod, da_obs, percentile=0.1)
        mae_1pctl = mean_absolute_error(da_mod, da_obs, percentile=0.01)
        rmse = root_mean_square_error(da_mod, da_obs)
        corr = spearman_correlation(da_mod, da_obs)
        binwidth = get_userdefined_binwidth(variable)  # Replace `None` with actual variable if needed
        pss, _, _, _ = perkins_skill_score(da_mod, da_obs, binwidth=binwidth)

        # Create a dictionary with the metrics and their values
        metrics_data = {
            "metric": [
                "Mean Bias",
                "Mean Absolute Error",
                "MAE at 90th Percentile",
                "MAE at 99th Percentile",
                "MAE at 10th Percentile",
                "MAE at 1st Percentile",
                "Root Mean Square Error",
                "Spearman Correlation",
                "Perkins Skill Score"
            ],
            member: [
                bias,
                mae,
                mae_90pctl,
                mae_99pctl,
                mae_10pctl,
                mae_1pctl,
                rmse,
                corr,
                pss
            ]
        }

        # Create a DataFrame
        if i == 0:
            df_metric = pd.DataFrame(metrics_data)
        else: 
            df_metric = pd.merge(df_metric, pd.DataFrame(metrics_data))

    # Save as CSV file
    if save_csv:
        # No user-specified path, use default
        if csv_path is None:
            csv_path = f"../output/metrics_{variable}_{da_obs.dataset}_{da_mod.dataset}.csv"

        df_metric.set_index('metric').to_csv(csv_path)

    return df_metric

def get_ranks_metrics(df: pd.DataFrame): 
    """
    Ranks the performance of different models across various metrics based on predefined ranking criteria.

    This function applies custom ranking rules to evaluate the performance of models across different metrics.
    The ranking is based on the following criteria:
    
    - 'Mean Bias' is ranked by its absolute value, with smaller values (closer to zero) ranked higher.
    - 'Spearman Correlation' and 'Perkins Skill Score' are ranked in descending order, meaning higher values (closer to 1) are better.
    - All other metrics are ranked in ascending order, where lower values are better.
    
    The input DataFrame `df` is expected to have the following structure:
    - The first column contains the metric names.
    - Each subsequent column contains the performance values of different models for each metric.
    
    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame where each row corresponds to a metric, the first column is the metric name, 
        and the subsequent columns contain performance values for different models.
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame where each value is replaced by its rank based on the ranking criteria for the corresponding metric.
        The rows are indexed by the metric names.
    
    """

    # Function to rank values
    def rank_values(row):

        if row['metric'] == 'Mean Bias':
            return row[1:].abs().rank(ascending=True)

        if row['metric'] in ['Spearman Correlation', 'Perkins Skill Score']:
            return row[1:].rank(ascending=False, method='min')
        else:
            return row[1:].rank(ascending=True, method='min')

    # Apply ranking
    df_ranked = df.apply(rank_values, axis=1).set_index(df['metric'])

    return df_ranked

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
