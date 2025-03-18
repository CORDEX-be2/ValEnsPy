import xarray as xr
import numpy as np
from scipy.stats import spearmanr
from datatree import DataTree
import pandas as pd
from functools import partial

from valenspy.processing import select_point
from valenspy.diagnostic.wrappers import acceptable_variables, required_variables

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

def annual_cycle(ds: xr.Dataset):
    """Calculate the annual cycle of the data. If lat and lon are present, the data is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the annual cycle of.

    Returns
    -------
    xr.Dataset
        The annual cycle of the data.
    """
    ds = _average_over_dims(ds, ["lat", "lon"])

    return ds.groupby("time.month").mean("time")


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

def time_series_trend(ds: xr.Dataset, window_size, min_periods: int = None, center: bool = True, **window_kwargs):
    """Calculate the trend of the time series data. If lat and lon are present, the data is averaged over the spatial dimensions lat and lon.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the trend of.
    window_size : int
        The size - in number of time steps - of the window to use for the rolling average.
    min_periods : int, optional
        The minimum number of periods required for a value to be considered valid, by default None
    center : bool, optional
        If True, the value is placed in the center of the window, by default True

    Returns
    -------
    xr.Dataset
        The trend of the data.
    """
    return _average_over_dims(ds, ["lat", "lon"]).rolling(time=window_size, min_periods=min_periods, center=center, **window_kwargs).mean()

def spatial_time_mean(ds: xr.Dataset):
    """Calculate the spatial mean of the data. If the time dimension is present, the data is averaged over the time dimension.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the spatial mean of.

    Returns
    -------
    xr.Dataset
        The spatial mean of the data.
    """
    return _average_over_dims(ds, "time")

@acceptable_variables(["tas", "tasmax", "tasmin"])
def urban_heat_island(ds: xr.Dataset, urban_coord: tuple, rural_coord: tuple, projection=None):
    """
    Calculate the urban heat island effect as the difference in temperature between an urban and rural area. 
    The grid-boxes closest to the urban and rural coordinates are selected and a difference is calculated between the two.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the urban heat island effect of.
    urban_coord : tuple
        The coordinates of the urban area in the format (lat, lon).
    rural_coord : tuple
        The coordinates of the rural area in the format (lat, lon).
    projection : str, optional
        The projection used to convert the urban and rural coordinates to the dataset's projection.
    
    Returns
    -------
    xr.Dataset
        The urban heat island effect as the difference in temperature between the urban and rural area.
    """
    urban = select_point(ds, lat_point=urban_coord[0], lon_point=urban_coord[1], projection=projection)
    rural = select_point(ds, lat_point=rural_coord[0], lon_point=rural_coord[1], projection=projection)

    return urban - rural

@acceptable_variables(["tas", "tasmax", "tasmin"])
def urban_heat_island_diurnal_cycle(ds: xr.Dataset, urban_coord: tuple, rural_coord: tuple, projection=None):
    """
    Calculate the diurnal cycle of the urban heat island effect as the difference in temperature between an urban and rural area.
    The grid-boxes closest to the urban and rural coordinates are selected and a difference is calculated between the two.

    Parameters
    ----------
    ds : xr.Dataset
        The data to calculate the urban heat island effect of.
    urban_coord : tuple
        The coordinates of the urban area in the format (lat, lon).
    rural_coord : tuple
        The coordinates of the rural area in the format (lat, lon).
    projection : str, optional
        The projection used to convert the urban and rural coordinates to the dataset's projection.

    Returns
    -------
    xr.Dataset
        The diurnal cycle of the urban heat island effect as the difference in temperature between the urban and rural area.
    """
    return diurnal_cycle(urban_heat_island(ds, urban_coord, rural_coord, projection=projection))

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

def calc_metrics_da(da_mod: xr.DataArray, da_obs: xr.DataArray, metrics=None, pss_binwidth=None):
    """
    Calculate statistical performance metrics for model data against observed data for a single variable.
    """
    if metrics is None:
        if not pss_binwidth:
            binwidth = get_userdefined_binwidth(da_mod.name)
        else:
            binwidth = pss_binwidth
        metrics = {
            "mean_bias": mean_bias,
            "mean_absolute_error": mean_absolute_error,
            "mae_90pctl": partial(mean_absolute_error, percentile=0.9),
            "mae_99pctl": partial(mean_absolute_error, percentile=0.99),
            "mae_10pctl": partial(mean_absolute_error, percentile=0.1),
            "mae_1pctl": partial(mean_absolute_error, percentile=0.01),
            "rmse": root_mean_square_error,
            "spearman_correlation": spearman_correlation,
            "PSS": partial(perkins_skill_score_value, binwidth=binwidth)
        }

    return {metric: metrics[metric](da_mod, da_obs) for metric in metrics}  
 
def calc_metrics_ds(ds_mod: xr.Dataset, ds_obs: xr.Dataset, metrics=None, pss_binwidth=None):
    """
    Calculate statistical performance metrics for model data against observed data for a dataset.
    """
    return {variable: calc_metrics_da(ds_mod[variable], ds_obs[variable], metrics, pss_binwidth=pss_binwidth) for variable in ds_mod.data_vars}

#####################################
# Ensemble2Ref diagnostic functions #
#####################################

def calc_metrics_dt(dt_mod: DataTree, da_obs: xr.Dataset, metrics=None, pss_binwidth=None):
    """
    Calculate statistical performance metrics for model data against observed data.

    This function computes various metrics between the model data stored in the DataTree 
    object `dt_mod` and the observed data `da_obs`.
    Default metrics include Mean Bias, Mean Absolute Error (MAE) at different percentiles,
    Root Mean Square Error (RMSE), Spearman Correlation, and Perkins Skill Score (PSS).
    
    Parameters:
    -----------
    dt_mod : DataTree
        A DataTree containing the model data for different members. The function loops 
        through each member to calculate the metrics.
    da_obs : xr.DataSet
        The observed data to compare against the model data.
    metrics : dict, optional
        A dictionary containing the names of the metrics to calculate and the corresponding functions. Default is the below specified metrics.
    pss_binwidth : float, optional
        The bin width to use for the Perkins Skill Score (PSS) calculation. If not provided, the optimal bin width is calculated.a
    
    Returns:
    --------
    df_metric : pd.DataFrame
        A DataFrame containing the calculated metrics and corresponding rank per metric for each member and variable in the data tree.

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

    data = {member.name: calc_metrics_ds(member.ds, da_obs, metrics=metrics, pss_binwidth=pss_binwidth) for member in dt_mod.leaves}
    #Create a pandas dataframe from a dictionary of dictionaries - each unique value of the most inner dictionary is a row
    df = pd.DataFrame.from_dict({(i,j): data[i][j]
                             for i in data.keys()
                             for j in data[i].keys()},
                            orient='index').reset_index().rename(columns={'level_0': 'member', 'level_1': 'variable'})
    df = df.melt(id_vars=["member", "variable"], var_name="metric", value_name="value")
    df = _add_ranks_metrics(df)
    return df


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

def _add_ranks_metrics(df: pd.DataFrame):
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
    
    def rank_per_metric(df):
        if df['metric'].iloc[0] == "mean_bias":
            df["rank"] = df["value"].abs().rank(ascending=True, method='min')
        elif df['metric'].iloc[0] in ['spearman_correlation', 'PSS']:
            df["rank"] = df["value"].rank(ascending=False, method='min')
        else:
            df["rank"] = df["value"].rank(ascending=True, method='min')
        return df

    df = df.groupby(["variable", "metric"]).apply(rank_per_metric).reset_index(drop=True)
    return df


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

def perkins_skill_score_value(da: xr.DataArray, ref: xr.DataArray, binwidth: float = None):
    return perkins_skill_score(da, ref, binwidth)[0]

