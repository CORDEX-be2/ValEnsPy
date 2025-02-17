# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
# %matplotlib inline

# %%
import xarray as xr
from datatree import DataTree
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

import cartopy.crs as ccrs
import valenspy as vp
from valenspy.diagnostic.visualizations import _add_features
from valenspy.diagnostic.functions import root_mean_square_error

git_dir = Path(os.popen("git rev-parse --show-toplevel").read().strip())

# %%
from valenspy.input.unit_converter import CORDEX_VARIABLES

# %% [markdown]
# User options

# %%
# Load data options
variables = ["tas"]  # "pr", "tasmin", "tasmax"
period = [1980, 2019]

# Plotting options - to be moved to a seperate "config" file
## Use the seaborn-v0_8-deep style?

mpl.rc('axes',edgecolor='grey')
mpl.rc('axes',labelcolor='dimgrey')
mpl.rc('axes',titlesize=14)
mpl.rc('axes',labelsize=10)
mpl.rc('axes', titlelocation="right")
mpl.rc('xtick',color='dimgrey')
mpl.rc('xtick',labelsize=10)
mpl.rc('ytick',color='dimgrey')
mpl.rc('ytick',labelsize=10)
mpl.rc('legend',fontsize='large')
mpl.rc('legend',frameon=False)
mpl.rc('text',color='dimgrey')
mpl.rc('figure', titlesize=16)


## Specify colors for the specific models
color_dict = {
    "/RCM/ERA5/ALARO1_SFX"          : "blue",
    "/RCM/ERA5/CCLM6-0-1-URB-ESG"   : "red",
    "/RCM/ERA5/MAR"                 : "green",
    "/obs/CLIMATE_GRID"             : "black"
}

d_cmap_diverging = { "tas": 'RdBu_r', "tasmax":'RdBu_r', "tasmin":'RdBu_r', "pr": 'BrBG' }

# Diagnostic switches

## Model2Self
do_AnnualCycle = True
do_TimeSeries = True
do_TimeSeriesUkkel = True
do_Trends = True
do_SpatialMean = True
## Model2Ref
do_SpatialBias = True


#Add iteration over the variables that we would like to be plotted.
#Add color maps per diagnostic
#Add option for season

# %% [markdown]
# STEP 1: Load the data

# %%
manager = vp.InputManager(machine="hortense")

# ALARO (Using the catalog and the variables userinput)

df_alaro = pd.read_csv("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/CORDEX_eval_scripts/catalog.csv")
df_alaro = df_alaro[df_alaro['frequency'] == 'day']
df_alaro = df_alaro[df_alaro['variable_id'].isin(variables)]
df_alaro

#Load all the paths in the df into one xarray dataset
ds_alaro = xr.open_mfdataset(df_alaro['path'].values, decode_coords='all', chunks="auto")
ds_alaro

## (Requires user adjustment)
# COSMO (Using the input manager - currently variables loaded manually)
experiment      = "CB2_CCLM_BEL28_ERA5_evaluation"
ds_cclm_tas     = manager.load_data("CCLM", ["tas"], freq="daily", path_identifiers=[experiment, "mean"])
# ds_cclm_tasmax  = manager.load_data("CCLM", ["tas"], freq="daily", path_identifiers=[experiment, "max"]).rename({'tas':'tasmax'})
# ds_cclm_tasmin  = manager.load_data("CCLM", ["tas"], freq="daily", path_identifiers=[experiment, "min"]).rename({'tas':'tasmin'})
# ds_cclm_pr      = manager.load_data("CCLM", ["pr"], freq="daily", path_identifiers=[experiment, "sum"])
# ds_cclm         = xr.merge([ds_cclm_tas, ds_cclm_pr, ds_cclm_tasmax, ds_cclm_tasmin])
#del ds_cclm_tas , ds_cclm_pr, ds_cclm_tasmax, ds_cclm_tasmin
ds_cclm = ds_cclm_tas

## (Requires user adjustment)
# MAR (Placeholder for MAR data - for plotting purposes)
ds_mar = ds_alaro

# Observational data
## CLIMATE_GRID (Regridded data)

ds_ref = manager.load_data("CLIMATE_GRID", variables, path_identifiers=["regridded"])

# %% [markdown]
# Create a DataTree object

# %%
## (Requires user adjustment)
data_dict = {
    "RCM/ERA5/ALARO1_SFX": ds_alaro,
    "RCM/ERA5/CCLM6-0-1-URB-ESG": ds_cclm,
    "RCM/ERA5/MAR": ds_mar,
    "obs/CLIMATE_GRID": ds_ref
}

dt = DataTree.from_dict(data_dict)

# %% [markdown]
# STEP 2: Preprocessing the data

# %%
#Regid to CLIMATE_GRID
dt["RCM"] = dt["RCM"].map_over_subtree(vp.remap_xesmf, dt.obs.CLIMATE_GRID.to_dataset(), method="conservative", regridding_kwargs={"keep_attrs": True})

#Select the time period from period[0] to period[1] (inclusive)
dt = dt.sel(time=slice(f"{period[0]}-01-01", f"{period[1]}-12-31"))

# %% [markdown]
# STEP 3: Diagnostics

# %%
from valenspy.diagnostic import AnnualCycle
from valenspy.diagnostic import TimeSeriesSpatialMean
from valenspy.diagnostic import TimeSeriesTrendSpatialMean
from valenspy.diagnostic import SpatialBias

# %% [markdown]
# Model2Self

# %% [markdown]
# Annual Cycle

# %%
if do_AnnualCycle:
    with ProgressBar():
        dt_annual_cycle = AnnualCycle(dt).compute()

# %%
if AnnualCycle:
    fig, ax = plt.subplots(figsize=(15, 5))
    AnnualCycle.plot_dt(dt_annual_cycle, var="tas", ax=ax, label="name", colors=color_dict)
    plt.legend()
    plt.savefig(git_dir / "CORDEX_eval_scripts/plots/Bel_mean_annual_cycle.png")

# %% [markdown]
# TimeSeries

# %%
if do_TimeSeries:
    heat_wave=["1997-06-01", "1997-08-31"]
    with ProgressBar():
        dt_time_series = dt.sel(time=slice(heat_wave[0], heat_wave[1]))
        dt_time_series = TimeSeriesSpatialMean(dt_time_series).compute()

# %%
if do_TimeSeries:
    fig, ax = plt.subplots(figsize=(15, 5))
    TimeSeriesSpatialMean.plot_dt(dt_time_series, var="tas", ax=ax, label="name", colors=color_dict)
    plt.legend()
    plt.savefig(git_dir / "CORDEX_eval_scripts/plots/Bel_mean_heat_wave.png")

# %% [markdown]
# TimeSeries (Ukkel)

# %%
if do_TimeSeriesUkkel:
    Ukkel = (4.37, 50.79)
    heat_wave=["1997-06-01", "1997-08-31"]
    with ProgressBar():
        dt_ukkel_hw = dt.map_over_subtree(vp.select_point, Ukkel[0], Ukkel[1]).sel(time=slice(heat_wave[0], heat_wave[1]))
        dt_time_series_uccle = TimeSeriesSpatialMean(dt_ukkel_hw).compute()

# %%
if do_TimeSeriesUkkel:
    fig, ax = plt.subplots(figsize=(15, 5))
    TimeSeriesSpatialMean.plot_dt(dt_time_series_uccle, var="tas", ax=ax, label="name", colors=color_dict)
    plt.legend()
    plt.savefig(git_dir / "CORDEX_eval_scripts/plots/Ukkel_heat_wave.png")

# %% [markdown]
# Trends

# %%
if do_Trends:
    Ukkel = (4.37, 50.79)
    dt_Ukkel = dt.map_over_subtree(vp.select_point, Ukkel[0], Ukkel[1])
    dt_Ukkel = dt_Ukkel - dt_Ukkel.sel(time=slice("1980-01-01", "1985-12-31")).mean("time")
    with ProgressBar():
        dt_trends = TimeSeriesTrendSpatialMean(dt_Ukkel, window_size=366*5).compute()

# %%
if do_Trends:
    fig, ax = plt.subplots(figsize=(15, 5))
    TimeSeriesTrendSpatialMean.plot_dt(dt_trends, var="tas", ax=ax, label="name", colors=color_dict)
    plt.legend()
    plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/Uccle_time_series_trend_normalized.png")

# %% [markdown]
# Model2Ref

# %% [markdown]
# SpatialBias

# %%
if do_SpatialBias:
    with ProgressBar():
        ds_ref = dt.obs.CLIMATE_GRID.to_dataset()
        dt_spatial_bias = SpatialBias(dt.RCM, ref=ds_ref).compute()

# %%
if do_SpatialBias:
    for var in variables:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": ccrs.PlateCarree()})
        axes = axes.flatten()
        cbar_kwargs={"shrink": 0.52}

        min = np.max([dt[var].values for dt in dt_spatial_bias.max().leaves])
        max = np.min([dt[var].values for dt in dt_spatial_bias.min().leaves])

        v_value = np.max([np.abs(min), np.abs(max)])

        for ax, dt_leave in zip(axes, dt_spatial_bias.leaves):
            SpatialBias.plot(dt_leave[var], ax=ax, cmap=d_cmap_diverging[var], cbar_kwargs=cbar_kwargs, vmin=-v_value, vmax=v_value)
            ax.set_title(dt_leave.name)
            #ax.set_title(' ')
            _add_features(ax, region='belgium')

        fig.suptitle(f'{CORDEX_VARIABLES[var]["long_name"]} bias compared to CLIMATE_GRID', y=0.8)
        fig.tight_layout()
        plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/spatialbias_{var}.png")

# %%
if SpatialBias:
#Plot option 2
    for var in variables:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": ccrs.PlateCarree()})
        axes = axes.flatten()
        cbar_kwargs={"shrink": 0.52}

        SpatialBias.plot_type = "facetted"

        SpatialBias.plot_dt(dt_spatial_bias, var=var, axes=axes, shared_cbar="abs", label="name", cbar_kwargs=cbar_kwargs, cmap=d_cmap_diverging[var])
        
        for ax in axes:
            _add_features(ax, region='belgium')

        fig.suptitle(f'{CORDEX_VARIABLES[var]["long_name"]} bias compared to CLIMATE_GRID', y=0.8)
        fig.tight_layout()
        plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/spatialbias_{var}.png")

# %%
  
## SpatialMean
if do_SpatialMean:
    d_cmap_sequential = {
    "tas": 'YlOrRd',
    "tasmax": 'YlOrRd',
    "tasmin": 'YlOrRd',
    "pr": 'YlGnBu' }
    with ProgressBar():
        dt_spatial_mean = dt.mean(dim="time")
        dt_spatial_mean = dt_spatial_mean.compute()

    #User input (size of the plots)
    fig, axes = plt.subplots(2,2,figsize=(9,5), subplot_kw={"projection": ccrs.PlateCarree()})
    axes = axes.flatten()

    # get overall min and max
    overall_min, overall_max = get_dt_overall_min_max(dt_timmean, variable)

    # loop over datasets to do plotting
    for i, dt_leave in enumerate(dt_timmean.leaves):
        ax = axes[i]
        ds = dt_leave.to_dataset()
        ds[variable].plot(ax=ax, cmap=d_cmap_sequential[variable], vmin=overall_min, vmax=overall_max)
        ax.set_title(dt_leave.name, loc='right')
        ax.set_title(' ')
        _add_features(ax, region='belgium')

    fig.suptitle(f'{ds[variable].long_name} ({ds[variable].units})', fontsize=16)

    fig.tight_layout()
    # Mean maps
    plt.savefig(f"/dodrio/scratch/users/vsc31332/rmi/vsc31332_inne/ValEnsPy/CORDEX_eval_scripts/plots/timmean_{variable}_{averaging_period}.png")


# %%
## SpatialMean
if SpatialMean:
    d_cmap_sequential = {
    "tas": 'YlOrRd',
    "tasmax": 'YlOrRd',
    "tasmin": 'YlOrRd',
    "pr": 'YlGnBu' }
    with ProgressBar():
        dt_spatial_mean = dt.mean(dim="time")
        dt_spatial_mean = dt_spatial_mean.compute()

    #User input (size of the plots)
    fig, axes = plt.subplots(2,2,figsize=(9,5), subplot_kw={"projection": ccrs.PlateCarree()})
    axes = axes.flatten()

    # get overall min and max
    overall_min, overall_max = get_dt_overall_min_max(dt_timmean, variable)

    # loop over datasets to do plotting
    for i, dt_leave in enumerate(dt_timmean.leaves):
        ax = axes[i]
        ds = dt_leave.to_dataset()
        ds[variable].plot(ax=ax, cmap=d_cmap_sequential[variable], vmin=overall_min, vmax=overall_max)
        ax.set_title(dt_leave.name, loc='right')
        ax.set_title(' ')
        _add_features(ax, region='belgium')

    fig.suptitle(f'{ds[variable].long_name} ({ds[variable].units})', fontsize=16)

    fig.tight_layout()
    # Mean maps
    plt.savefig(f"/dodrio/scratch/users/vsc31332/rmi/vsc31332_inne/ValEnsPy/CORDEX_eval_scripts/plots/timmean_{variable}_{averaging_period}.png")

