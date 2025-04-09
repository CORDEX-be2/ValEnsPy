import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
from valenspy._utilities._regions import region_bounds
from valenspy.diagnostic.functions import perkins_skill_score
from valenspy.diagnostic.plot_utils import default_plot_kwargs, _augment_kwargs

import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np

# make sure xarray passes the attributes when doing operations - change default for this
xr.set_options(keep_attrs=True)

###################################
# Model2Self diagnostic functions #
###################################


def plot_diurnal_cycle(da: xr.DataArray, **kwargs):
    """Plot the daily cycle of the data.
    
    Parameters
    ----------
    da : xr.DataArray
        The data array to plot the daily cycle of.
    **kwargs : dict
        Additional keyword arguments to pass to the xarray DataArray plot method.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes with the plotted daily cycle.
    """
    da.plot(**kwargs)

    ax = _get_gca(**kwargs)
    
    ax.set_xlabel('Hour of day')
    ax.set_xticks(range(0, 24, 3))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: '{:02d}:00'.format(int(x))))

    return ax

def plot_annual_cycle(da: xr.DataArray, **kwargs):
    """Plot the annual cycle of the data.

    Parameters
    ----------
    da : xr.DataArray
        The data array to plot the annual cycle of.
    **kwargs : dict
        Additional keyword arguments to pass to the xarray DataArray plot method.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes with the plotted annual cycle.
    """
    da.plot(**kwargs)

    ax = _get_gca(**kwargs)

    ax.set_xlabel('Month')

    return ax


def plot_time_series(da: xr.DataArray, **kwargs):
    """Plot a time series of the data.
    
    Parameters
    ----------
    da : xr.DataArray
        The data array to plot the time series of.
    **kwargs : dict
        Additional keyword arguments to pass to the xarray DataArray plot method.
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes with the plotted time series.
    """
    da.plot(**kwargs)

    ax = _get_gca(**kwargs)

    return ax

@default_plot_kwargs({
    'subplot_kws': {'projection': ccrs.PlateCarree()}
    })
def plot_map(da: xr.DataArray, **kwargs):
    """
    Plots a simple map of a 2D xarray DataArray.

    Default plot settings:

    - subplot_kws: {'projection': ccrs.PlateCarree()}
    - cbar_kwargs: {'label': "{{da.attrs.get('long_name', 'Data')}} ({{da.name}})"}
    

    Parameters:
    -----------
    da : xr.DataArray
        The 2D xarray DataArray to plot. It should have latitude and longitude dimensions.
    title : str, optional
        The title for the plot. If not provided, a default title based on the DataArray's
        long_name attribute will be set.
    **kwargs :
        Additional keyword arguments to pass to the xarray DataArray plot method.

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The matplotlib Axes with the plot.
    """
    kwargs = _augment_kwargs({"cbar_kwargs": {"label": f"{da.attrs.get('long_name', 'Data')} ({da.name})"}}, **kwargs)

    da.plot(**kwargs)

    ax = _get_gca(**kwargs)

    return ax


##################################
# Model2Ref diagnostic visuals   #
##################################

def create_custom_cmap(hex_color1: str, hex_color2: str, num_colors: int):
    """
    Create a custom colormap that transitions between two given hex colors and generates a specified number of colors.

    Parameters
    ----------
    hex_color1 : str
        The starting color of the colormap in hex format (e.g., '#FF0000' for red).
    hex_color2 : str
        The ending color of the colormap in hex format (e.g., '#0000FF' for blue).
    num_colors : int
        The number of colors to generate in the colormap, including both the start and end colors.

    Returns
    -------
    cmap : matplotlib.colors.LinearSegmentedColormap
        A colormap that can be used in plotting functions to visualize data with a color gradient from hex_color1 to hex_color2.
    colors : numpy.ndarray
        An array of the RGB values for each of the colors in the generated colormap.
    
    Example
    -------
    >>> cmap, colors = create_custom_cmap("#FF0000", "#0000FF", 10)
    >>> plt.imshow([colors], aspect='auto')
    >>> plt.show()

    """
    # Convert hex colors to RGB
    rgb_color1 = mcolors.hex2color(hex_color1)
    rgb_color2 = mcolors.hex2color(hex_color2)

    # Create colormap
    cmap = LinearSegmentedColormap.from_list("custom_cmap", [rgb_color1, rgb_color2], N=num_colors)


    return cmap


def plot_maps_mod_ref_diff(
    da_mod: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    region=None,
    **kwargs,
):
    """
    Plots side-by-side comparison maps for model data, reference data, and their difference.

    Parameters
    ----------
    da_mod : xarray.DataArray
        The model data to be plotted. Expected to be a 2D array with latitude and longitude dimensions.
    da_ref : xarray.DataArray
        The reference data to be plotted. Expected to be a 2D array with latitude and longitude dimensions.
    da_diff : xarray.DataArray
        The difference between model and reference data (model - reference) to be plotted.
    region : str, optional
        A string specifying the region for determining the plotting extent. Defaults to None.
    **kwargs : dict, optional
        Additional keyword arguments for customizing the plots:

            - `vmin` (float): Minimum value for the color scale of the model and reference plots.
            - `vmax` (float): Maximum value for the color scale of the model and reference plots.
            - `vmin_bias` (float): Minimum value for the color scale of the bias plot.
            - `vmax_bias` (float): Maximum value for the color scale of the bias plot.

            
    Returns
    -------
    list of matplotlib.axes._subplots.AxesSubplot
        A list containing the three axes objects for the model, reference, and bias plots.

    Notes
    -----

    - All warnings are suppressed during the execution of this function.
    - Each subplot uses the PlateCarree projection and includes borders and coastlines.
    - The color scale for the model and reference plots is determined by the combined minimum and maximum values
    across both datasets unless explicitly specified via `vmin` and `vmax`.
    - The bias plot uses a diverging color map ('coolwarm') centered around zero, with its range determined
    by the maximum absolute value of the difference unless overridden by `vmin_bias` and `vmax_bias`.
    - Titles for the subplots are positioned to the right and are derived from the 'dataset' attribute of the
    respective data arrays, if available. Otherwise, default titles ("Model", "Reference", "Model - Reference")
    are used.
    - The figure title is set to the 'long_name' attribute of the reference data, followed by its variable name.

    """

    # Turn off all warnings
    warnings.filterwarnings("ignore")

    fig, axes = plt.subplots(
        1, 3, figsize=(14, 3), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    axes = axes.flatten()

    # find plotting min and max
    cbar_label = f"{da_ref.attrs['long_name']} ({da_ref.attrs['units']})"
    cbar_kwargs = {"label": cbar_label}

    if "vmin" in kwargs:
        vmin = kwargs.pop("vmin")
    else:  # plotting boundaries
        vmin = float(min(da_mod.min().values, da_ref.min().values))

    if "vmax" in kwargs:
        vmax = kwargs.pop("vmax")
    else:  # plotting boundaries
        vmax = float(max(da_mod.max().values, da_ref.max().values))

    # titles - use the dataset attribute if available
    if "dataset" in da_mod.attrs:
        mod_title = da_mod.attrs["dataset"]
    else:
        mod_title = "Model"

    if "dataset" in da_ref.attrs:
        ref_title = da_ref.attrs["dataset"]
    else:
        ref_title = "Reference"

    # mod
    ax = axes[0]
    da_mod.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
    ax.set_title("")
    ax.set_title(mod_title, loc="right")
    _add_features(ax, region=region)

    # ref
    ax = axes[1]
    da_ref.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
    ax.set_title("")
    ax.set_title(ref_title, loc="right")
    _add_features(ax, region=region)

    # bias
    diff_bound = float(max(abs(da_diff.min().values), abs(da_diff.max().values)))

    if "vmin_bias" in kwargs:
        vmin = kwargs.pop("vmin_bias")
    else:  # plotting boundaries
        vmin = -diff_bound

    if "vmax_bias" in kwargs:
        vmax = kwargs.pop("vmax_bias")
    else:  # plotting boundaries
        vmax = diff_bound

    ax = axes[2]
    da_diff.plot(
        ax=ax,
        cmap="coolwarm",
        vmax=vmax,
        vmin=-diff_bound,
        cbar_kwargs=cbar_kwargs,
    )

    ax.set_title("")
    ax.set_title(f"{mod_title} - {ref_title}", loc="right")
    _add_features(ax, region=region)

    fig.suptitle(f"{da_ref.attrs['long_name']} ({da_ref.name})", y=1)
    fig.tight_layout()

    return axes


def plot_time_series_mod_ref(
    da_mod: xr.DataArray, da_ref: xr.DataArray, ax=None, title: str = None, **kwargs
):
    """
    Plot time series for both model and reference datasets on the same axes.

    Parameters
    ----------
    da_mod : xarray.DataArray
        The DataArray containing the model time series data.
    da_ref : xarray.DataArray
        The DataArray containing the reference time series data.
    ax : matplotlib.axes.Axes, optional
        The axes on which to plot the time series. If None, a new figure and axes are created.
    title : str, optional
        The title for the plot. If None, a default title based on `da_mod` attributes is used.
    **kwargs : dict
        Additional keyword arguments passed to `xarray.DataArray.plot`.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted time series.
    """
    if ax is None:
        fig, ax = plt.subplots()

    # Plot the reference data array on the axes
    da_ref.plot(ax=ax, label=da_ref.attrs.get("dataset", "Reference"), color="k")

    # Plot the model data array on the same axes with some transparency
    da_mod.plot(ax=ax, label=da_mod.attrs.get("dataset", "Model"), alpha=0.5, **kwargs)

    # Add a legend without a frame
    ax.legend(frameon=False)

    # Set the title, either the provided one or based on the model data attributes
    if title is None:
        ax.set_title(f"{da_mod.attrs.get('long_name', 'Data')} ({da_mod.name})")
    else:
        ax.set_title(title)

    return ax


def plot_points_on_map(d_point_coords: dict, ax=None, region=None):
    """
    Plot geographic points on a map using Cartopy, with optional region highlighting.

    Parameters
    ----------
    d_point_coords : dict
        A dictionary where keys are point identifiers (e.g., station names or IDs) and values are tuples of
        longitude and latitude coordinates (e.g., {'Point1': (lon1, lat1), 'Point2': (lon2, lat2)}).
    ax : matplotlib.axes.Axes, optional
        The axes on which to plot the points. If None, a new figure and axes with a PlateCarree projection are created.
    region : str or None, optional
        The region to highlight on the map. This could be a predefined region name (e.g., 'belgium')
        or None if no specific region is needed.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted points and the map features.

    Example
    -------
    >>> d_point_coords = {'Point1': (4.3517, 50.8503), 'Point2': (5.5413, 50.6326)}
    >>> plot_points_on_map(d_point_coords, region="belgium")
    """
    # Create a figure and set the projection to PlateCarree
    if ax is None:
        fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

    # Plot each point and add a label
    for point_id, (lon, lat) in d_point_coords.items():
        ax.plot(
            lon,
            lat,
            marker="o",
            color="red",
            markersize=5,
            transform=ccrs.PlateCarree(),
        )
        ax.text(lon + 0.1, lat - 0.1, point_id, transform=ccrs.PlateCarree())

    # Add coastline and country borders and region selection if region is provided
    _add_features(ax, region=region)

    ax.set_title("Location of points", loc="right")

    return ax


def visualize_perkins_skill_score(da_mod: xr.DataArray, da_obs: xr.DataArray, binwidth: float = None):
    """
    Visualize the Perkins Skill Score (PSS) by plotting the normalized histograms 
    of the model and reference data, and display the PSS score and bin width used.
    For testing bin_widths
    
    Parameters
    ----------
    da_mod : xr.DataArray
        The model data to compare.
    da_obs : xr.DataArray
        The reference data to compare against.
    binwidth : float, optional
        The width of each bin for the histogram. If None, an optimal bin width 
        should be calculated within the function (default is None).
    
    Returns
    -------
    None
        This function does not return any value. It displays a plot with the 
        normalized histograms and Perkins Skill Score.
    
    Notes
    -----
    The function calculates the Perkins Skill Score using the provided or default
    bin width, and plots the normalized histograms of the model and reference data.
    The plot also includes annotations for the Perkins Skill Score and bin width used.
    """
    # Calculate Perkins Skill Score and histograms
    pss_score, freq_m, freq_r, binwidth = perkins_skill_score(da_mod, da_obs, binwidth=binwidth)

    # Create the plot
    fig, ax = plt.subplots()

    # Plot the histograms
    ax.plot(freq_m, label="model")
    ax.plot(freq_r, label="ref", color="k")
    ax.set_title('Normalized histograms for calculating Perkins Skill Score', loc='right')
    ax.set_xlabel('bins')
    ax.set_ylabel('frequency')
    ax.legend(frameon=False, loc='upper right')

    # Annotate the plot with PSS score and bin width
    ax.text(0.05, 0.9, f"Perkins skill score: {pss_score:.3f}", transform=ax.transAxes)
    ax.text(0.05, 0.85, f"Used binwidth: {binwidth:.2f}", transform=ax.transAxes)

    # Adjust layout
    fig.tight_layout()
    plt.show()

def plot_metric_ranking(df_metric, ax=None, plot_colorbar=True, hex_color1 = None, hex_color2=None, **kwargs):
    """
    Plots a heatmap of the ranking of metrics for different model members.

    This function takes a DataFrame of metrics, calculates the rankings of these metrics 
    for each model member, and creates a heatmap representing the ranks. The plot can 
    optionally include a colorbar to represent the ranking levels. If no axis is provided, 
    a new figure and axis are created for the plot.

    Parameters:
    -----------
    df_metric : pd.DataFrame
        A DataFrame containing the calculated metrics for different model members. Each column represents a model member, and each row represents a metric.
    ax : matplotlib.axes.Axes, optional
        A pre-existing axis to plot the heatmap. If None (default), a new figure and axis 
        are created.
    plot_colorbar : bool, optional
        If True (default), a colorbar is added to the plot to represent the rank levels. 
        If False, the heatmap is plotted without a colorbar.
    hex_color1 : str
        The starting color of the colormap in hex format (e.g., '#FF0000' for red).
    hex_color2 : str
        The ending color of the colormap in hex format (e.g., '#0000FF' for blue).

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axis object containing the heatmap plot.

    Notes:
    ------
    - The color map has the 'summer' palette as default and is resampled to the number of model members.
    - A customized color map can be included or determined as an interpolation between two colorcodes (hex codes)
    - Rankings are normalized based on the number of model members.
    - The function supports colorbar ticks to represent custom rank labels, which are added 
      only if `plot_colorbar=True`.
    """

    df_p = df_metric.pivot(index=["metric"], columns="member", values="rank")
    df_p = df_p[df_p.sum().sort_values().index]

    num_levels = df_metric["member"].nunique()
    if "cmap" not in kwargs:
        if hex_color1 and hex_color2:
            cmap = create_custom_cmap(hex_color1=hex_color1, hex_color2=hex_color2, num_colors=num_levels)
        else:
            cmap = plt.get_cmap('summer', num_levels)
    else:
        cmap = plt.get_cmap(kwargs.pop("cmap"), num_levels)
    
    boundaries = np.arange(1, num_levels + 2, 1)
    norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

    if ax is None: 
        fig, ax = plt.subplots()

    if "title" in kwargs:
        ax.set_title(kwargs.pop("title"), loc="right")

    heatmap = sns.heatmap(df_p, ax=ax, cbar=plot_colorbar, cmap=cmap, norm=norm, **kwargs)
    ax.set_ylabel(' ')
    ax.set_xlabel('Members')
    
    if plot_colorbar:
        colorbar = heatmap.collections[0].colorbar
        colorbar.set_ticks(np.arange(1, num_levels + 1) + .5)  # Set the ticks you want
        colorbar.set_ticklabels(range(1, num_levels + 1))  # Set the custom labels for the ticks
    
    return ax

##################################
# Helper functions               #
##################################

def _get_gca(**kwargs):
    """
    Get the current axes as the axis passed in the keyword arguments or the current axis.
    """
    if "ax" in kwargs:
        return kwargs["ax"]
    else:
        return plt.gca()

# Define a function to add borders, coastlines to the axes
def _add_features(ax, region=None):
    """
    Adds geographical features to a given cartopy GeoAxes.

    Parameters:
    ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The GeoAxes to which features are to be added.

    Features Added:
    - Borders: Adds country borders with a dotted linestyle.
    - Coastlines: Adds coastlines with a specified linewidth and color.
    - extent: if region is given, cut out the plotting extent based on the lat and lon bounds given.

    Notes:
    - The function can be extended to set the extent of the plot by uncommenting and modifying the
      set_extent line to include appropriate longitude and latitude bounds.

    Example:
    >>> fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    >>> add_features(ax)
    >>> plt.show()
    """

    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color="k")

    if not region is None:
        lon_bounds = region_bounds[region]["lon_bounds"]
        lat_bounds = region_bounds[region]["lat_bounds"]
        ax.set_extent(
            [lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]],
            crs=ccrs.PlateCarree(),
        )
