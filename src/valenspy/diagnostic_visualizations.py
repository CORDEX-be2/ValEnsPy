import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
from valenspy._regions import region_bounds

# make sure xarray passes the attributes when doing operations - change default for this
xr.set_options(keep_attrs=True)

###################################
# Model2Self diagnostic functions #
###################################


def plot_diurnal_cycle(data: xr.DataArray, ax, **kwargs):
    """Plot the daily cycle of the data."""
    data.plot(ax=ax, **kwargs)
    ax.set_title("Daily Cycle")
    return ax


def plot_time_series(da: xr.DataArray, ax=None, **kwargs):
    """
    Plot a time series from an xarray DataArray.

    Parameters
    ----------
    da : xarray.DataArray
        The DataArray containing the time series data to plot.
    ax : matplotlib.axes.Axes, optional
        The axes on which to plot the time series. If None, a new figure and axes are created.
    **kwargs : dict
        Additional keyword arguments passed to `xarray.DataArray.plot`.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted time series.
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    # Plot the data array on the provided or newly created axes
    da.plot(ax=ax, **kwargs)
    
    # Set the title based on the 'long_name' attribute
    ax.set_title(da.attrs.get('long_name', ''), loc='left')
    ax.set_title(' ', loc='center')

    return ax


def plot_map(da: xr.DataArray, ax=None, title=None, region=None, **kwargs):
    
    """
    Plots a simple map of a 2D xarray DataArray.

    This function creates a map plot for a given 2D xarray DataArray, optionally using 
    a provided matplotlib Axes, which needs to have a projection. It automatically sets the colorbar label based on the 
    DataArray's attributes and adds features like coastlines and country borders.

    Parameters:
    -----------
    da : xr.DataArray
        The 2D xarray DataArray to plot. It should have latitude and longitude dimensions.
    ax : matplotlib.axes.Axes, optional
        The matplotlib Axes on which to plot the map. If not provided, a new Axes with a 
        PlateCarree projection will be created.
    title : str, optional
        The title for the plot. If not provided, a default title based on the DataArray's 
        long_name attribute will be set.
    **kwargs : 
        Additional keyword arguments to pass to the xarray DataArray plot method.

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The matplotlib Axes with the plot.

    Example:
    --------
    >>> import xarray as xr
    >>> import matplotlib.pyplot as plt
    >>> import cartopy.crs as ccrs
    >>> fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    >>> ax = plot_map(da, ax=ax)
    """
    
    # If ax is not provided, create a new one
    if ax is None:
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Set colorbar label
    cbar_kwargs = {'label': f"{da.attrs.get('long_name', 'Data')} ({da.attrs.get('units', '')})"}

    # Plot the data array with the specified colorbar axis
    da.plot(ax=ax, cbar_kwargs=cbar_kwargs, **kwargs)

    if title is None: 
        # Set the title
        ax.set_title(f"{da.attrs.get('long_name', 'Data')} ({da.name})")

    # Add coastline and country borders and region selection if region is provided
    _add_features(ax, region=region)

    return ax

##################################
# Model2Ref diagnostic visuals   #
##################################



def plot_spatial_bias(da: xr.DataArray, ax=None, region = None, **kwargs):
    """Plot the spatial bias of the data compared to the reference."""

    # if no ax element is passed, create one
    if ax is None: 
        fig , ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    

    title = f"Mean bias of {da.long_name}"
    cmap = "coolwarm"
    plot_map(da, ax=ax, title=title, cmap = cmap)
    ax.set_title(title)
    # Add coastline and country borders and region selection if region is provided
    _add_features(ax, region=region)

    return ax


def plot_maps_mod_ref_diff(da_mod: xr.DataArray,  da_ref: xr.DataArray,  da_diff: xr.DataArray, region=None, return_fig=False): 

  """
  Plots comparison maps for model data, reference data, and their difference.

  Parameters:
  da_mod (xarray.DataArray): The model data to be plotted: 2D lat
  da_ref (xarray.DataArray): The reference data to be plotted.
  da_diff (xarray.DataArray): The difference (model - reference) to be plotted.
  region (str)              : string of the region to determine the plotting extend.
  return_fig (boolean)      : determines whether the figure object is returned, default False

  Returns (optionally):
  matplotlib.figure.Figure: The figure object containing the three subplots.

  Notes:
  - This function suppresses all warnings.
  - Each subplot is created with the PlateCarree projection and includes borders and coastlines.
  - The color scale for the model and reference plots is determined by the minimum and maximum values 
    across both datasets.
  - The bias plot uses a diverging color map ('coolwarm') centered around zero.
  - The color bar label is derived from the 'long_name' and 'units' attributes of the reference data.
  - The titles of the subplots are positioned to the right.
  - The figure title is set to the 'long_name' attribute of the reference data.
  
  """


  # Turn off all warnings
  warnings.filterwarnings('ignore')

  fig, axes = plt.subplots(1, 3, figsize=(14, 3), subplot_kw={'projection': ccrs.PlateCarree()})
  axes = axes.flatten()

  # find plotting min and max
  cbar_label = f"{da_ref.attrs['long_name']} ({da_ref.attrs['units']})"
  cbar_kwargs = {'label': cbar_label}


  # plotting bounds
  vmin = float(min(da_mod.min().values, da_ref.min().values))
  vmax = float(max(da_mod.max().values, da_ref.max().values))


  # titles - use the dataset attribute if available
  if 'dataset' in da_mod.attrs: 
      mod_title = da_mod.attrs['dataset']
  else: 
      mod_title = 'Model'

  if 'dataset' in da_ref.attrs: 
      ref_title = da_ref.attrs['dataset']
  else: 
      ref_title = 'Reference'

  # mod
  ax = axes[0]
  da_mod.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
  ax.set_title('')
  ax.set_title(mod_title, loc='right')
  _add_features(ax, region=region)

  # ref
  ax = axes[1]
  da_ref.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
  ax.set_title('')
  ax.set_title(ref_title, loc='right')
  _add_features(ax, region=region)

  # bias
  ax = axes[2]
  diff_bound = float(max(abs(da_diff.min().values), abs(da_diff.max().values)))
  da_diff.plot(ax=ax, cmap = 'coolwarm', vmax = diff_bound, vmin = - diff_bound, cbar_kwargs=cbar_kwargs)
  ax.set_title('')
  ax.set_title(f"{mod_title} - {ref_title}", loc='right')
  _add_features(ax, region=region)

  fig.suptitle(f"{da_ref.attrs['long_name']} ({da_ref.name})", y=1);
  fig.tight_layout()

  if return_fig: 
    return fig


def plot_time_series_mod_ref(da_mod: xr.DataArray, da_ref: xr.DataArray, ax=None, title: str = None, **kwargs):
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
    da_ref.plot(ax=ax, label=da_ref.attrs.get("dataset", "Reference"), color='k')
    
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
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Plot each point and add a label
    for point_id, (lon, lat) in d_point_coords.items():
        ax.plot(lon, lat, marker='o', color='red', markersize=5, transform=ccrs.PlateCarree())
        ax.text(lon + 0.1, lat - 0.1, point_id, transform=ccrs.PlateCarree())

    # Add coastline and country borders and region selection if region is provided
    _add_features(ax, region=region)
    
    ax.set_title('Location of points', loc='right')

    return ax

##################################
# Helper functions               #
##################################



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


    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE, linewidth = 0.5, color = 'k')

    if not region is None: 
        lon_bounds = region_bounds[region]['lon_bounds']
        lat_bounds = region_bounds[region]['lat_bounds']
        ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())

