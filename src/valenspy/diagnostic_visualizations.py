import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings


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


##################################
# Model2Ref diagnostic visuals   #
##################################


def plot_spatial_bias(data: xr.DataArray, ax, **kwargs):
    """Plot the spatial bias of the data compared to the reference."""
    data.plot(
        ax=ax, cmap="coolwarm", cbar_kwargs={"label": "Temperature bias (K)"}, **kwargs
    )
    ax.set_title("Spatial Bias")

    return ax


def plot_maps_mod_ref_diff(
    da_mod: xr.DataArray, da_ref: xr.DataArray, da_diff: xr.DataArray
):
    """
    Plots comparison maps for model data, reference data, and their difference.

    Parameters:
    da_mod (xarray.DataArray): The model data to be plotted: 2D lat
    da_ref (xarray.DataArray): The reference data to be plotted.
    da_diff (xarray.DataArray): The difference (model - reference) to be plotted.

    Returns:
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
    warnings.filterwarnings("ignore")

    fig, axes = plt.subplots(
        1, 3, figsize=(14, 3), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    axes = axes.flatten()

    # find plotting min and max
    cbar_label = f"{da_ref.attrs['long_name']} ({da_ref.attrs['units']})"
    cbar_kwargs = {"label": cbar_label}

    # plotting bounds
    vmin = float(min(da_mod.min().values, da_ref.min().values))
    vmax = float(max(da_mod.max().values, da_ref.max().values))

    # titles
    mod_title = da_mod.attrs["dataset"]
    ref_title = da_ref.attrs["dataset"]

    # mod
    ax = axes[0]
    da_mod.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
    ax.set_title("")
    ax.set_title(mod_title, loc="right")
    _add_features(ax)

    # ref
    ax = axes[1]
    da_ref.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
    ax.set_title("")
    ax.set_title(ref_title, loc="right")
    _add_features(ax)

    # bias
    ax = axes[2]
    diff_bound = float(max(abs(da_diff.min().values), abs(da_diff.max().values)))
    da_diff.plot(
        ax=ax,
        cmap="coolwarm",
        vmax=diff_bound,
        vmin=-diff_bound,
        cbar_kwargs=cbar_kwargs,
    )
    ax.set_title("")
    ax.set_title(f"{mod_title} - {ref_title}", loc="right")
    _add_features(ax)

    fig.suptitle(f"{da_ref.attrs['long_name']} ({da_ref.name})", y=1)
    fig.tight_layout()
    plt.show()

    # fig.savefig(f"./plots/{region}_{experiment}_{ref_dataset}_timmean_bias.png")
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

    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color="k")
    # ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
