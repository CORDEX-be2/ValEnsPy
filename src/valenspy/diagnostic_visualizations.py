import xarray as xr
import matplotlib.pyplot as plt

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


# Standard matplotlib settings
# Setting domains for the data - coastlines, cartopy, achtegronden.
