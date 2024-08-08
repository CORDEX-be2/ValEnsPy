import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings

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


def plot_time_series(data: xr.DataArray, ax, **kwargs):
    data.plot(ax=ax, **kwargs)
    ax.set_title("Time Series")
    return ax

def plot_map(da: xr.DataArray, ax=None, title=None, **kwargs):
    
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

    # Add coastline and country borders
    _add_features(ax)

    return ax

##################################
# Model2Ref diagnostic visuals   #
##################################



def plot_spatial_bias(da: xr.DataArray, ax=False, **kwargs):
    """Plot the spatial bias of the data compared to the reference."""

    # if no ax element is passed, create one
    if not ax: 
        fig , ax= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    

    title = f"Mean bias of {da.long_name}"
    cmap = "coolwarm"
    plot_map(da, ax=ax, title=title, cmap = cmap)

    return ax


def plot_maps_mod_ref_diff(da_mod: xr.DataArray,  da_ref: xr.DataArray,  da_diff: xr.DataArray, return_fig=False): 

  """
  Plots comparison maps for model data, reference data, and their difference.

  Parameters:
  da_mod (xarray.DataArray): The model data to be plotted: 2D lat
  da_ref (xarray.DataArray): The reference data to be plotted.
  da_diff (xarray.DataArray): The difference (model - reference) to be plotted.
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
  _add_features(ax)

  # ref
  ax = axes[1]
  da_ref.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
  ax.set_title('')
  ax.set_title(ref_title, loc='right')
  _add_features(ax)

  # bias
  ax = axes[2]
  diff_bound = float(max(abs(da_diff.min().values), abs(da_diff.max().values)))
  da_diff.plot(ax=ax, cmap = 'coolwarm', vmax = diff_bound, vmin = - diff_bound, cbar_kwargs=cbar_kwargs)
  ax.set_title('')
  ax.set_title(f"{mod_title} - {ref_title}", loc='right')
  _add_features(ax)

  fig.suptitle(f"{da_ref.attrs['long_name']} ({da_ref.name})", y=1);
  fig.tight_layout()

  if return_fig: 
    return fig


##################################
# Helper functions               #
##################################



# Define a function to add borders, coastlines to the axes
def _add_features(ax):

    """
    Adds geographical features to a given cartopy GeoAxes.

    Parameters:
    ax (cartopy.mpl.geoaxes.GeoAxesSubplot): The GeoAxes to which features are to be added.

    Features Added:
    - Borders: Adds country borders with a dotted linestyle.
    - Coastlines: Adds coastlines with a specified linewidth and color.

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
    #ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())

