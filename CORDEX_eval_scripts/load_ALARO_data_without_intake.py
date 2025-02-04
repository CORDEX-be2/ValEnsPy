import xarray as xr
import valenspy as vp
from datatree import DataTree
from pathlib import Path
import pandas as pd
from dask.diagnostics import ProgressBar

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

dt.to_netcdf("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/notebooks/intermediate_data/ALARO_interm_20250204.nc")