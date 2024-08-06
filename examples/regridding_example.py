import valenspy as vp
from valenspy._utilities import load_yml
import xarray as xr
from pathlib import Path
import cdo

# Illustration of simple usage of the vp.preprocessing_tasks regridding functionality.
# This script will regrid the CCLM data to the ERA5 grid using the bilinear method.
#################
# Load the data
#################

# ERA5 data
manager = vp.InputManager(machine="hortense")
ds_ref = manager.load_data(
    "ERA5",
    "tas",
    period=[1995, 1995],
    freq="daily",
    region="europe",
    path_identifiers=["-daily-"],
)

# CCLM data to be put into a input converter and manager
## CCLM input manager will produce the paths
dataset_PATHS = load_yml("dataset_PATHS")
model_directory = dataset_PATHS["hortense"]["CCLM"]
mod_LOOKUP = load_yml(f"CCLM_lookup")

# get CCLM variable corresponding to the requested variable using its look-up table
mod_var = mod_LOOKUP["tas"]["mod_name"]

# define the path
directory = Path(model_directory + "EUR11_CO_TA_GC_TSO" + "/" + mod_var + "/")

# define the CCLM files for the corresponding variable
mod_files = list(directory.glob(mod_var + "_daymean.nc"))  # Se
ds_mod = xr.open_mfdataset(mod_files, combine="by_coords", chunks="auto")


#################
# Regridding
#################

gridfile = manager._get_file_paths(
    "ERA5",
    "tas",
    period=[1995, 1995],
    freq="daily",
    region="europe",
    path_identifiers=["-daily-"],
)[0]

from valenspy.preprocessing_tasks.regrid import remap_cdo

ds_mod_regrid = remap_cdo(gridfile, ds_mod, remap_method="bil")

ds_mod_regrid = ds_mod_regrid.rename({"T_2M": "tas"})
#################
# Plotting
#################
# Take the average over the time dimension of both datasets, plot this and then plot the difference between the two.
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ds_ref["tas"].mean("time").plot(ax=ax[0])
ds_mod_regrid["tas"].mean("time").plot(ax=ax[1])
(ds_ref["tas"] - ds_mod_regrid["tas"]).mean("time").plot(ax=ax[2])
plt.show()
fig.savefig("ERA5_CCLM_bias.png")
