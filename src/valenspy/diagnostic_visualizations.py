import xarray as xr
import matplotlib.pyplot as plt

##################################
# Model2Ref diagnostic visuals   #
##################################

def plot_spatial_bias(data: xr.DataArray, ax, **kwargs):
    data.plot(
        ax=ax, cmap="coolwarm", cbar_kwargs={"label": "Temperature bias (K)"}, **kwargs
    )
    ax.set_title("Spatial Bias")

    return ax

# Standard matplotlib settings
# Setting domains for the data - coastlines, cartopy, achtegronden.
