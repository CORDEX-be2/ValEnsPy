import xarray as xr
import matplotlib.pyplot as plt
##################################
# Model2Ref diagnostic visuals   #
##################################

def plot_spatial_bias(data: xr.DataArray, **kwargs):
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    data.plot(ax=ax, cmap="coolwarm", cbar_kwargs={"label": "Temperature bias (K)"}, **kwargs)
    ax.set_title("Spatial Bias")

    return fig

#Standard matplotlib settings
#Setting domains for the data - coastlines, cartopy, achtegronden.
