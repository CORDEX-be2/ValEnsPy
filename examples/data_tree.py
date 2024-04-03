import sys, os
from pathlib import Path
import numpy

# lib_folder = Path(__file__).resolve().parents[1].joinpath("src")

# test_data_folder = Path(__file__).resolve().parents[1].joinpath("tests").joinpath("data")
#Get the git root directory

test_data_folder = r"C:\Users\kvandela\Projects\ValEnsPy\tests\data"

#Import xarray datatree
from datatree import DataTree

import xarray as xr
# Create a DataTree object
ds=xr.open_dataset(Path(test_data_folder).joinpath("tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc"))

ds=xr.open_mfdataset(test_data_folder + "/*historical*.nc", combine='by_coords', chunks='auto')
ds2=xr.open_mfdataset(test_data_folder + "/*ssp245*.nc", combine='by_coords', chunks='auto')

dt = DataTree.from_dict({"historical": ds, "ssp245": ds2})

avg=dt.mean("time")
avg

#Plot the average temperature for the historical and ssp245 scenario
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
avg.ssp245.tas.plot(ax=ax[0], label="ssp245")
avg.historical.tas.plot(ax=ax[1], label="historical")
#Plot the difference between the two scenarios
(avg.ssp245.tas-avg.historical.tas).plot(ax=ax[2], label="ssp245 - historical")
plt.show()

