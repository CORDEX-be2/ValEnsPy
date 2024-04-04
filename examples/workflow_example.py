from pathlib import Path

import valenspy as vp
import xarray as xr
from datatree import DataTree

#First make an inputconverter object - assume that all the data is already in the correct format (which it will likely not be)
from valenspy.inputconverter_functions import _non_convertor

ic = vp.InputConverter(_non_convertor)

EOBS_data_dir = Path("/dodrio/scratch/projects/2022_200/project_input/External/observations/EOBS/0.1deg/")

EOBS_obs_files = list(EOBS_data_dir.glob("*.nc")) #Select all the netCDF files in the directory
EOBS_obs_files = ic.convert_input(EOBS_obs_files) #Convert the input to the correct format

EOBS_ds = xr.open_mfdataset(EOBS_obs_files[0:2], combine='by_coords', chunks='auto')

#Now make an ensemble member object
EC_Earth3_dir = Path("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/tests/data")
EC_Earth3_hist_files = list(EC_Earth3_dir.glob("*historical*.nc")) #Select all the netCDF files in the directory

EC_Earth3_hist_files = ic.convert_input(EC_Earth3_hist_files) #Convert the input to the correct format

EC_Earth3_ds = xr.open_mfdataset(EC_Earth3_hist_files, combine='by_coords', chunks='auto')

#Now make an "future" ensemble member object
EC_Earth3_ssp_files = list(EC_Earth3_dir.glob("*ssp245*.nc")) #Select all the netCDF files in the directory
EC_Earth3_ssp_files = ic.convert_input(EC_Earth3_ssp_files) #Convert the input to the correct format
EC_Earth3_ssp_ds = xr.open_mfdataset(EC_Earth3_ssp_files, combine='by_coords', chunks='auto')

#Now make a datatree object
dt = DataTree.from_dict({"obs/EOBS": EOBS_ds, "ensembles/cmip6/EC_Earth3/hist": EC_Earth3_ds, "ensembles/cmip6/EC_Earth3/ssp245": EC_Earth3_ssp_ds})

#Apply some postprocessing operations on the datatree
pp = vp.Preprocessor()
pp.add_preprocessing_task(vp.preprocessing_tasks.Regrid(dt.obs.EOBS.ds, name="to_obs", description="Regrid the CMIP6 data to the EOBS grid"))
pp.add_preprocessing_task(vp.preprocessing_tasks.Area(dt.obs.EOBS.ds, name="to_obs", description="Select the common area of the EOBS grid"))

#Apply the preprocessing
dt.ensembles = pp.apply_preprocessing(dt.ensembles)

#Note that all ensemble (model data) objects in the datatree have been regridded to the EOBS grid
dt

#Apply a diagnostic operation
#First make a diagnostic object
from valenspy.diagnostic_functions import spatial_bias
from valenspy.diagnostic_visualizations import plot_spatial_bias

diag = vp.Model2Obs(spatial_bias, plot_spatial_bias, name="spatial_bias", description="Calculate the time averaged spatial bias between the model and the observations")

#Apply the diagnostic
ds = diag.apply(dt.ensembles.cmip6.EC_Earth3.hist, dt.obs.EOBS)

