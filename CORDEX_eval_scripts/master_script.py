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

# %% [markdown]
# # CORDEX.be II evaluation 
#

# %%
# %matplotlib inline

import xarray as xr
from datatree import DataTree, map_over_subtree
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import xclim

import cartopy.crs as ccrs
import valenspy as vp
from valenspy.diagnostic.visualizations import _add_features
from valenspy.diagnostic.functions import root_mean_square_error
from valenspy.input.unit_converter import CORDEX_VARIABLES


git_dir = Path(os.popen("git rev-parse --show-toplevel").read().strip())

# %% [markdown]
# ## User options

# %%
# Load data options
variables = ["tas", "pr", "tasmin", "tasmax"]
period = [1985, 2019]

#Unit handeling - choose the units for the output
unit_dict = {
    "tas": "C",
    "tasmax": "C",
    "tasmin": "C",
    "pr": "mm/day"
}

# Plotting options - see matplotlibrc file

## Specify colors for the specific models
color_dict = {
    "/RCM/ERA5/ALARO1_SFX"          : "tab:blue",
    "/RCM/ERA5/CCLM6-0-1-URB-ESG"   : "tab:orange",
    "/RCM/ERA5/MAR"                 : "tab:green",
    "/obs/CLIMATE_GRID"             : "black"
}

d_cmap_diverging = { "tas": 'RdBu_r', "tasmax":'RdBu_r', "tasmin":'RdBu_r', "pr": 'BrBG', "tnn":'RdBu_r', "txx":'RdBu_r', "rx1day": 'BrBG'}
d_cmap_sequential = {"tas": 'YlOrRd', "tasmax": 'YlOrRd', "tasmin": 'YlOrRd', "pr": 'YlGnBu', "tnn":'YlOrRd', "txx":'YlOrRd', "rx1day": 'YlGnBu'}

#Diagnostic options
## Model2Self
do_AnnualCycle = {
    "compute": True,
    "variables" : variables,
}
do_TimeSeries = {
    "compute": True,
    "variables" : variables,
    "periods": [["1997-06-01", "1997-08-31"], ["2003-07-20", "2003-08-20"]]
}
do_TimeSeriesUkkel = {
    "compute": True,
    "variables" : variables,
    "periods": [["1997-06-01", "1997-08-31"]]
}
do_Trends = {
    "compute": True,
    "variables" : variables
}
do_SpatialMean = {
    "compute": True,
    "variables" : variables,
    "reference": "/obs/CLIMATE_GRID",
    "seasons" : ["All", "DJF", "MAM", "JJA", "SON"]
}

## Model2Ref
do_SpatialBias = {
    "compute": True,
    "variables" : variables,
    "reference": "/obs/CLIMATE_GRID",
    "seasons" : ["All" , "DJF", "MAM", "JJA", "SON"]
}

#Add xclim variables

# %% [markdown]
# ## STEP 1: Load the data

# %%
manager = vp.InputManager(machine="hortense")

# ALARO (Using the catalog and the variables userinput)

df_alaro = pd.read_csv(git_dir/"CORDEX_eval_scripts/catalog.csv")
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
ds_cclm_tasmax  = manager.load_data("CCLM", ["tas"], freq="daily", path_identifiers=[experiment, "max"]).rename({'tas':'tasmax'})
ds_cclm_tasmin  = manager.load_data("CCLM", ["tas"], freq="daily", path_identifiers=[experiment, "min"]).rename({'tas':'tasmin'})
ds_cclm_pr      = manager.load_data("CCLM", ["pr"], freq="daily", path_identifiers=[experiment, "sum"])
ds_cclm         = xr.merge([ds_cclm_tas, ds_cclm_pr, ds_cclm_tasmax, ds_cclm_tasmin])

# adjust the time from 11:30 to 12:00 to make xclim work. 
new_time = ds_cclm.tas.time.astype('datetime64[D]') + np.timedelta64(12, 'h')
ds_cclm = ds_cclm.assign_coords(time=new_time)

del ds_cclm_tas , ds_cclm_pr, ds_cclm_tasmax, ds_cclm_tasmin


## (Requires user adjustment)
# MAR (Placeholder for MAR data - for plotting purposes)

ds_mar = xr.open_mfdataset("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc44757_Nicolas/remap_mar/MAR-ERA5/*.nc")
from valenspy.input import INPUT_CONVERTORS
ds_mar = INPUT_CONVERTORS["MAR"].convert_input(ds_mar)
#Keep only the variables of interest
ds_mar = ds_mar[variables]

# Observational data
## CLIMATE_GRID (Regridded data)

ds_ref = manager.load_data("CLIMATE_GRID", variables, path_identifiers=["regridded"])


# %%

# %% [markdown]
# ### Create a DataTree object

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
# ## STEP 2: Preprocessing the data

# %% [markdown]
# ### CCLM xclim to be fixed - freq attribute?

# %%
#To be added in valenspy (including a non map_over_subtree version) - which is basically a wrapper over xclim.units.convert_units_to
@map_over_subtree
def convert_units_to(ds, var, target_unit, context="infer"):
    if var in ds:
        ds = ds.copy()
        ds[var] = xclim.units.convert_units_to(ds[var], target_unit, context=context)
    return ds


# %%
#Move to valenspy when happy with the implementation
#How to deal with indicators that require multiple variables?s
@map_over_subtree
def xclim_indicator(ds, indicator, vars, **kwargs):
    if isinstance(vars, str):
        return indicator(ds[vars], **kwargs).to_dataset()
    elif isinstance(vars, list): #Order is important!
        data_arrays = [ds[var] for var in vars]
        return indicator(*data_arrays, **kwargs).to_dataset()



# %%
#Regid to CLIMATE_GRID
dt["RCM"] = dt["RCM"].map_over_subtree(vp.remap_xesmf, dt.obs.CLIMATE_GRID.to_dataset(), method="conservative", regridding_kwargs={"keep_attrs": True})

#Select the time period from period[0] to period[1] (inclusive)
dt = dt.sel(time=slice(f"{period[0]}-01-01", f"{period[1]}-12-31"))

#Unit conversion to the desired units
for var, unit in unit_dict.items():
    dt = convert_units_to(dt, var, unit)

dt = dt.compute()

# %%
#Working but requries CCLM data fix
# Location is important! Before or after regridding? Before or after unit conversion? Best before!
dt_tnn = xclim_indicator(dt, xclim.indicators.cf.tnn, vars="tas", freq="YS")
dt_txx = xclim_indicator(dt, xclim.indicators.cf.txx, vars="tasmax", freq="YS")
dt_max_1day_precipitation_amount = xclim_indicator(dt, xclim.atmos.max_1day_precipitation_amount, vars="pr", freq="YS")

xclim_dt_dict = {
     "tnn": dt_tnn,
     "txx": dt_txx,
     "rx1day": dt_max_1day_precipitation_amount
 }

# %% [markdown]
# ## STEP 3: Diagnostics

# %%
from valenspy.diagnostic import AnnualCycle
from valenspy.diagnostic import TimeSeriesSpatialMean
from valenspy.diagnostic import TimeSeriesTrendSpatialMean
from valenspy.diagnostic import SpatialTimeMean
from valenspy.diagnostic import SpatialBias

# %% [markdown]
# ### Model2Self

# %% [markdown]
# #### Annual Cycle

# %%
if do_AnnualCycle["compute"]:
    with ProgressBar():
        print("Computing annual cycle")
        dt_annual_cycle = AnnualCycle(dt).compute()

    for var in do_AnnualCycle["variables"]:
        fig, ax = plt.subplots(figsize=(10, 6))
        AnnualCycle.plot_dt(dt_annual_cycle, var=var, ax=ax, label="name", colors=color_dict)
        plt.title(f"Annual cycle of {CORDEX_VARIABLES[var]['long_name']}")
        plt.legend()
        plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_bel_mean_annual_cycle.png",bbox_inches='tight')
        plt.close(fig)

# %% [markdown]
# #### TimeSeries

# %%
if do_TimeSeries["compute"]:
    for period in do_TimeSeries["periods"]:
        if period == "All":
            dt_time_series = dt
        else:
            dt_time_series = dt.sel(time=slice(period[0], period[1]))
        with ProgressBar():
            print(f"Computing time series for period {period}")
            dt_time_series = TimeSeriesSpatialMean(dt_time_series).compute()

        for var in do_TimeSeries["variables"]:
            fig, ax = plt.subplots(figsize=(15, 5))
            TimeSeriesSpatialMean.plot_dt(dt_time_series, var=var, ax=ax, label="name", colors=color_dict)
            plt.title(f"Time series of {CORDEX_VARIABLES[var]['long_name']}")
            plt.legend()
            period_filename = "_".join(period)
            plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_time_series_bel_mean_{period_filename}.png",bbox_inches='tight')
            plt.close(fig)

# %% [markdown]
# ##### Xclim variables and derivatives

# %%
for var, dt_xclim in xclim_dt_dict.items():
     with ProgressBar():
         print(f"Computing TimeSeriesSpatialMean for {var}")
         dt_time_series = TimeSeriesSpatialMean(dt_xclim).compute()

     fig, ax = plt.subplots(figsize=(15, 5))
     TimeSeriesSpatialMean.plot_dt(dt_time_series, var=var, ax=ax, label="name", colors=color_dict)
     plt.title(f"Time series of {var}")
     plt.legend()
     plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_time_series_bel_mean.png",bbox_inches='tight')
     plt.close(fig)

# %% [markdown]
# #### TimeSeries (Ukkel)

# %%

if do_TimeSeriesUkkel["compute"]:
    Ukkel = (4.37, 50.79)
    dt_ukkel_hw = dt.map_over_subtree(vp.select_point, lat_point=Ukkel[1], lon_point=Ukkel[0])
    for period in do_TimeSeries["periods"]:
        if period == "All":
            dt_time_series = dt_ukkel_hw
        else:
            dt_time_series = dt_ukkel_hw.sel(time=slice(period[0], period[1]))
        
        with ProgressBar():
            print(f"Computing time series for {period}")
            dt_time_series_uccle = TimeSeriesSpatialMean(dt_time_series).compute()

        for var in do_TimeSeries["variables"]:
            fig, ax = plt.subplots(figsize=(15, 5))
            TimeSeriesSpatialMean.plot_dt(dt_time_series_uccle, var=var, ax=ax, label="name", colors=color_dict)
            plt.title(f"Time series of {CORDEX_VARIABLES[var]['long_name']} in Ukkel")
            plt.legend()
            plt.tight_layout()
            period_filename = "_".join(period)
            plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_time_series_ukkel_{period_filename}.png",bbox_inches='tight')
            plt.close(fig)

# %% [markdown]
# #### Trends

# %%
if do_Trends["compute"]:
    years=5
    window_size = 366*years # 5 years
    Ukkel = (4.37, 50.79)
    dt_Ukkel = dt.map_over_subtree(vp.select_point, Ukkel[0], Ukkel[1])
    dt_Ukkel = dt_Ukkel - dt_Ukkel.sel(time=slice("1980-01-01", "1985-12-31")).mean("time")
    with ProgressBar():
        dt_trends = TimeSeriesTrendSpatialMean(dt_Ukkel, window_size=window_size).compute()

    for var in do_Trends["variables"]:
        fig, ax = plt.subplots(figsize=(15, 5))
        TimeSeriesTrendSpatialMean.plot_dt(dt_trends, var=var, ax=ax, label="name", colors=color_dict)
        plt.legend()
        plt.title(f"Trend of {CORDEX_VARIABLES[var]['long_name']} in Ukkel")
        plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_trend_ukkel_window_{years}y.png",bbox_inches='tight')
        plt.close(fig)

# %% [markdown]
# ##### Xclim variables and derivatives

# %%
for var, dt_xclim in xclim_dt_dict.items():
    years=5
    window_size = years # 5 years
    Ukkel = (4.37, 50.79)
    dt_Ukkel = dt_xclim.map_over_subtree(vp.select_point, Ukkel[0], Ukkel[1])
    dt_Ukkel = dt_Ukkel - dt_Ukkel.sel(time=slice("1980-01-01", "1985-12-31")).mean("time")
    with ProgressBar():
        dt_trends = TimeSeriesTrendSpatialMean(dt_Ukkel, window_size=window_size).compute()

    fig, ax = plt.subplots(figsize=(15, 5))
    TimeSeriesTrendSpatialMean.plot_dt(dt_trends, var=var, ax=ax, label="name", colors=color_dict)
    plt.legend()
    plt.title(f"Trend of {var} in Ukkel")
    plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_trend_ukkel_window_{years}y.png",bbox_inches='tight')
    plt.close(fig)

# %% [markdown]
# #### Spatial Mean

# %%
if do_SpatialMean["compute"]:
    #Add other processing steps here (e.g. period selection)
    for season in do_SpatialMean["seasons"]:
        #Select season
        if season != "All":
            dt_s = dt.map_over_subtree(lambda x: x.sel(time=x.time.dt.season == season))
        else:
            dt_s = dt
        #Compute
        with ProgressBar():
            print(f"Computing spatial mean for {season}")
            dt_spatial_mean = SpatialTimeMean(dt_s).compute()

        #Plot
        for var in do_SpatialMean["variables"]:
            fig, axes = plt.subplots(2, 2, figsize=(9, 5), subplot_kw={"projection": ccrs.PlateCarree()})
            axes = axes.flatten()
            
            SpatialTimeMean.plot_type = "facetted"

            SpatialTimeMean.plot_dt(dt_spatial_mean, var=var, axes=axes,shared_cbar="min_max", label="name", cmap=d_cmap_sequential[var])

            for ax in axes:
                _add_features(ax, region='belgium')
            
            fig.suptitle(f'{CORDEX_VARIABLES[var]["long_name"]} spatial mean (Season: {season})',  y=1.01)
            fig.tight_layout()
            plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_spatial_mean_{season}.png",bbox_inches='tight')
            plt.close(fig)

# %% [markdown]
# ##### Xclim variables and derivatives

# %%
for var, dt_xclim in xclim_dt_dict.items():
    with ProgressBar():
        print(f"Computing spatial mean for {var}")
        dt_spatial_mean = SpatialTimeMean(dt_xclim).compute()

    fig, axes = plt.subplots(2, 2, figsize=(9, 5), subplot_kw={"projection": ccrs.PlateCarree()})
    axes = axes.flatten()

    SpatialTimeMean.plot_type = "facetted"

    SpatialTimeMean.plot_dt(dt_spatial_mean, var=var, axes=axes, shared_cbar="min_max", label="name", cmap=d_cmap_sequential[var])

    for ax in axes:
        _add_features(ax, region='belgium')

    fig.suptitle(f'{var} spatial mean',  y=1.01)
    fig.tight_layout()
    plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_spatial_mean_{season}.png",bbox_inches='tight')


# %% [markdown]
# ### Model2Ref

# %% [markdown]
# #### SpatialBias

# %%
if do_SpatialBias["compute"]:
    for season in do_SpatialBias["seasons"]:
        if season != "All":
            dt_s = dt.map_over_subtree(lambda x: x.sel(time=x.time.dt.season == season))
        else:
            dt_s = dt

        #Compute    
        with ProgressBar():
            print(f"Computing spatial bias for season {season}")
            ds_ref = dt_s[do_SpatialBias["reference"]].to_dataset()
            dt_spatial_bias = SpatialBias(dt_s.RCM, ref=ds_ref).compute()

        #Plot
        for var in do_SpatialBias["variables"]:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": ccrs.PlateCarree()})
            axes = axes.flatten()
            cbar_kwargs={"shrink": 0.52}

            SpatialBias.plot_dt(dt_spatial_bias, var=var, axes=axes, shared_cbar="abs", label="name", cbar_kwargs=cbar_kwargs, cmap=d_cmap_diverging[var])
            
            for ax in axes:
                _add_features(ax, region='belgium')

            fig.suptitle(f'{CORDEX_VARIABLES[var]["long_name"]} bias compared to CLIMATE_GRID (Season: {season})', y=0.8)
            fig.tight_layout()
            plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_spatialbias_{season}.png",bbox_inches='tight')
            plt.close(fig)


# %% [markdown]
# ##### Xclim variables and derivatives

# %%
for var, dt_xclim in xclim_dt_dict.items():
    with ProgressBar():
        print(f"Computing spatial bias for {var}")
        ds_ref = dt_xclim[do_SpatialBias["reference"]].to_dataset()
        dt_spatial_bias = SpatialBias(dt_xclim.RCM, ref=ds_ref).compute()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": ccrs.PlateCarree()})
    axes = axes.flatten()
    cbar_kwargs={"shrink": 0.52}

    SpatialBias.plot_dt(dt_spatial_bias, var=var, axes=axes, shared_cbar="abs", label="name", cbar_kwargs=cbar_kwargs, cmap=d_cmap_diverging[var])

    for ax in axes:
        _add_features(ax, region='belgium')

    fig.suptitle(f'{var} bias compared to CLIMATE_GRID', y=0.8)
    fig.tight_layout()
    plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/{var}_spatialbias.png",bbox_inches='tight')
    plt.close(fig)
