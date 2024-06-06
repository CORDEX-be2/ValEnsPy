#!/usr/bin/env python3

"""
Read in CLIMATE_GRID in .Rdata, plot on grid and save as netCDF

Reading in .RData files after retrieval from Oracle using script from Bert

"""

# %%

import rdata
import xarray as xr

variables = ["EVAPOTRANS_REF", "SUN_INT", "SUN_DURATION", "PRECIP_DURATION", "WIND_PEAK_SPEED",
	         "PRECIP_1H_MAX", "EVAPOTRANS_REF", "TEMP_MAX","HUMIDITY_RELATIVE","TEMP_MIN", "TEMP_AVG",
	         "WIND_SPEED", "PRESSURE", "SHORT_WAVE_FROM_SKY", "SUN_INT_HORIZ", "PRECIP_QUANTITY"]


variables = ["EVAPOTRANS_REF"]

start_yr = 1950
end_yr = 2023
data_dir = '/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc31332_inne/CLIMATE_GRID/'

for variable in variables: 

    filename = 'climate_atlas_'+str(variable)+'_CLIMATE_GRID_'+str(start_yr)+'_'+str(end_yr)+'_raw.RData'

# read in data as dict and extract pandas df
d_data = rdata.read_rda(data_dir+filename)

df = d_data['all.data']



# %% Create xarray data array

da = xr.DataArray(
    data=temperature,
    dims=["lat", "lon", "time"],
    coords=dict(
        lon=(["x", "y"], lon),
        lat=(["x", "y"], lat),
        time=time,
        reference_time=reference_time,
    ),
    attrs=dict(
        description="Ambient temperature.",
        units="degC",
    ),
)


filename_out = str(variable)'_CLIMATE_GRID_'+str(start_yr)+'_'+str(end_yr)+'_daily.nc"
