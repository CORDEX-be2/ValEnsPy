import xarray as xr
import valenspy as vp
from datatree import DataTree
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
import os

#Get the git directory using the Path object using an os command
git_dir = Path(os.popen("git rev-parse --show-toplevel").read().strip())


#User options
variables = ["tas", "pr"]
period = [1980,2019]
target_grid="/dodrio/scratch/projects/2022_200/external/climate_grid/TEMP_AVG_CLIMATE_GRID_1954_2023_daily.nc"

############################################
# STEP 1: Loading the data

# start input manager
manager = vp.InputManager(machine="hortense")


#MODEL data
# Load the ALARO data

df_alaro = pd.read_csv("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/CORDEX_eval_scripts/catalog.csv")
df_alaro = df_alaro[df_alaro['frequency'] == 'day']
df_alaro = df_alaro[df_alaro['variable_id'].isin(variables)]
df_alaro

#Load all the paths in the df into one xarray dataset
ds_alaro = xr.open_mfdataset(df_alaro['path'].values, decode_coords='all', chunks="auto")
ds_alaro

# Load the COSMO data
## Issue - variarble names are hardcoded due to daily statistic, necessary to find correct path with input manager
experiment = "CB2_CCLM_EUR11_ERA5_evaluation"
ds_cclm_tas = manager.load_data("CCLM", ["tas"], freq="daily", path_identifiers=[experiment, "mean"])
ds_cclm_pr  = manager.load_data("CCLM", ["pr"], freq="daily", path_identifiers=[experiment, "sum"])
ds_cclm = xr.merge([ds_cclm_tas, ds_cclm_pr])
del ds_cclm_tas, ds_cclm_pr


# Load the MAR data
## palceholder for MAR data -for plotting purposes
ds_mar = ds_alaro

#OBSERVATIONAL data

# Load CLIMATE_GRID
## Issue - currently can't load the ungridded CLIMATE_GRID data as there is no unique identifier for the data 

ds_ref = manager.load_data("CLIMATE_GRID", variables, path_identifiers=["regridded"])

# Create a DataTree object
data_dict = {
    "RCM/ERA5/ALARO1_SFX": ds_alaro,
    "RCM/ERA5/CCLM6-0-1-URB-ESG": ds_cclm,
    "RCM/ERA5/MAR": ds_mar,
    "obs/CLIMATE_GRID": ds_ref
}

dt = DataTree.from_dict(data_dict)

############################################
# STEP 2: Preprocessing the data

## Regrid (currently to CLIMATE_GRID)
dt["RCM"] = dt["RCM"].map_over_subtree(vp.remap_xesmf, dt.obs.CLIMATE_GRID.to_dataset(), method="conservative", regridding_kwargs={"keep_attrs": True})

## Select the time period from period[0] to period[1] (inclusive)
dt = dt.sel(time=slice(f"{period[0]}-01-01", f"{period[1]}-12-31"))

############################################
# STEP 3: Diagnostics
#Model2Self

#AnnualCycle
from valenspy.diagnostic import AnnualCycle

with ProgressBar():
    dt_annual_cycle = dt.map_over_subtree(AnnualCycle)
    dt_annual_cycle = dt_annual_cycle.compute()

fig, ax = plt.subplots(figsize=(15, 5))

for leaf in dt_annual_cycle.leaves:
    AnnualCycle.plot(leaf["tas"], ax=ax, label=leaf.name)

plt.legend()
plt.savefig(git_dir / "CORDEX_eval_scripts/plots/Bel_mean_annual_cycle.png")

#TimeSeries
from valenspy.diagnostic import TimeSeriesSpatialMean

heat_wave=["1997-06-01", "1997-08-31"]

dt_heat = dt.sel(time=slice(heat_wave[0], heat_wave[1]))

with ProgressBar():
    dt_time_series = dt_heat.map_over_subtree(TimeSeriesSpatialMean)
    dt_time_series = dt_time_series.compute()

fig, ax = plt.subplots(figsize=(15, 5))

for leaf in dt_time_series.leaves:
    TimeSeriesSpatialMean.plot(leaf["tas"], ax=ax, label=leaf.name)

plt.legend()
plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/Bel_mean_time_series_{heat_wave[0]}_{heat_wave[1]}.png")

#TimeSeries - Ukkel
Ukkel = (4.37, 50.79)

#Note this is beter done in the native grid (here we use the regridded grid)
dt_heat_uccle = dt_heat.map_over_subtree(vp.select_point, Ukkel)

with ProgressBar():
    dt_time_series_uccle = dt_heat_uccle.map_over_subtree(TimeSeriesSpatialMean)
    dt_time_series_uccle = dt_time_series_uccle.compute()

fig, ax = plt.subplots(figsize=(15, 5))
for leaf in dt_time_series_uccle.leaves:
    TimeSeriesSpatialMean.plot(leaf["tas"], ax=ax, label=leaf.name)

plt.legend()
plt.savefig(git_dir / f"CORDEX_eval_scripts/plots/Uccle_time_series_{heat_wave[0]}_{heat_wave[1]}.png")

#Trend figures ??
ds = dt.obs.CLIMATE_GRID.to_dataset()
ds_m = ds.mean(dim=["lat","lon"])
ds_t = ds_m.rolling(time=365*10).mean()

fig, ax = plt.subplots(figsize=(15, 5))
ds_t.tas.plot()
plt.savefig(git_dir / "test.png")