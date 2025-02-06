import xarray as xr
import valenspy as vp
from datatree import DataTree
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

#User options
variables = ["tas", "pr"]
period = [1980,2002]
target_grid="/dodrio/scratch/projects/2022_200/external/climate_grid/TEMP_AVG_CLIMATE_GRID_1954_2023_daily.nc"
############################################

#MODEL data
# Load the ALARO data
df = pd.read_csv("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/CORDEX_eval_scripts/catalog.csv")
df = df[df['frequency'] == 'day']
df = df[df['variable_id'].isin(variables)]
df

#Load all the paths in the df into one xarray dataset
ds = xr.open_mfdataset(df['path'].values, decode_coords='all', chunks="auto")
ds

# Load the COSMO data

# Load the MAR data

#Load the observational data
## Issue - currently can't load the ungridded CLIMATE_GRID data as there is no unique identifier for the data 
manager = vp.InputManager(machine="hortense")
ds_ref = manager.load_data("CLIMATE_GRID", variables, path_identifiers=["regridded"])

# Create a DataTree object
data_dict = {
    "RCM/ERA5/ALARO1_SFX": ds,
    "obs/CLIMATE_GRID": ds_ref
}

dt = DataTree.from_dict(data_dict)

#Some preprocessing steps
## Regrid (currently to CLIMATE_GRID)
dt["RCM"] = dt["RCM"].map_over_subtree(vp.remap_xesmf, dt.obs.CLIMATE_GRID.to_dataset(), method="bilinear", regridding_kwargs={"keep_attrs": True})

## Select the time period from 1980 to 2002 (inclusive)
dt = dt.sel(time=slice(f"{period[0]}-01-01", f"{period[1]}-12-31"))

#Apply diagnostics
#Model2Self
diag = vp.diagnostic.DiurnalCycle
with ProgressBar():
    ds_alaro = diag.apply(dt["RCM/ERA5/ALARO1_SFX"].to_dataset())
    ds_alaro = ds_alaro.compute()
    ds_obs = diag.apply(dt["obs/CLIMATE_GRID"].to_dataset())
    ds_obs = ds_obs.compute()

fig, ax = plt.subplots(1, 2, figsize=(15, 5))
diag.plot(ds_alaro["tas"], ax=ax[0], label="ALARO1_SFX")
diag.plot(ds_obs["tas"], ax=ax[1], label="CLIMATE_GRID")
#Add a legend
plt.legend()
plt.show()
plt.savefig("diurnal_cycle.png")

##Model2Ref
diag = vp.diagnostic.SpatialBias
with ProgressBar():
    ds_spbias = diag.apply(dt["RCM/ERA5/ALARO1_SFX"], dt["obs/CLIMATE_GRID"])
    ds_spbias = ds_spbias.compute()
    
diag.plot(ds_spbias)