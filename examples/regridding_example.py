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
ds_era5 = manager.load_data(
    "ERA5",
    "tas, pr",
    period=[1995, 1995],
    freq="daily",
    region="europe",
    path_identifiers=["-daily-"],
)

# CCLM data to be put into a input converter and manager
## CCLM input manager will produce the paths
ds_eobs = manager.load_data("EOBS",["tas", "pr"], path_identifiers = ["0.1deg",  "mean"])

# select 1995
ds_eobs = ds_eobs.sel(time=ds_eobs.time.dt.year.isin(1995))

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

ds_eobs_regrid = remap_cdo(gridfile, ds_eobs, remap_method="con")


#################
# Plotting
#################
# Take the average over the time dimension of both datasets, plot this and then plot the difference between the two.
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ds_era5["tas"].mean("time").plot(ax=ax[0])
ds_eobs_regrid["tas"].mean("time").plot(ax=ax[1])
(ds_era5["tas"].mean("time")-ds_eobs_regrid["tas"].mean("time")).plot(ax=ax[2])
fig.tight_layout()