#!/usr/bin/env python3

"""
Regrid CLIMATE_GRID to regular lon/lat grid 

Script to regrid CLIMATE_GRID from LAMBERT into own defined LAT/LON grid
! includes a user-defined regular lat-lon grid, as close as possible to the native CLIMATE_GRID
Uses CDO regridding. 

I. Vanderkelen, June 2024
"""


# load modules 
import xarray as xr
import os
from pathlib import Path
import numpy as np


# ---------------------------------------
# STEP 1. Define grid and write grid

# define number of latitude and longitude grid cells
# regrid to lat lon grid
# lat lon resolution, close to 5 km in Belgium - define grid based on this. 

gridname = 'latlon_be_5km'

xinc = resolution_lon
yinc = resolution_lat
nlats=45
nlons=60
xfirst=2.4
yfirst=49.5

# save grid info as text
# Create the content to write to the file
content = f"""gridtype=lonlat
xsize={nlons}
ysize={nlats}
xfirst={xfirst}
xinc={xinc}
yfirst={yfirst}
yinc={yinc}
"""

# Write the content to a file
with open(gridname+".txt", "w") as file:
    file.write(content)

# ---------------------------------------
# STEP 2. Do remapping


# User settings
dataset = "CLIMATE_GRID"
directory = Path('/mnt/HDS_CLIMATE/CLIMATE/CLIMATE_GRID/')

remap_method = "remapcon" # cdo method

variables = ["EVAPOTRANS_REF", "SUN_INT", "SUN_DURATION", "PRECIP_DURATION", "WIND_PEAK_SPEED", "PRECIP_1H_MAX", "EVAPOTRANS_REF", "TEMP_MAX","HUMIDITY_RELATIVE", "TEMP_AVG", "WIND_SPEED", "PRESSURE", "SHORT_WAVE_FROM_SKY", "SUN_INT_HORIZ", "PRECIP_QUANTITY", "TEMP_MIN"]

# create regrid folder 
directory_grid = str(directory)+'/regridded/'+gridname+'/'

if not os.path.exists(directory_grid):
    os.makedirs(directory_grid) 

for variable in variables: 
    files = list(directory.glob("*"+variable+"*.nc")) # Select all the netCDF files for the variable in the directory

    for file in files:  

        # if grid already exists, don't do remapping. 
        if not gridname in str(file): 
            filename_out = str(file).replace('.nc','_'+gridname+'.nc')

            print('regridding '+str(file))
            # do remapping using CDO
            os.system('cdo '+remap_method+','+gridname+".txt "+ str(file) + " "+filename_out)

            # move regridded file to regridded directory
            os.system('mv '+filename_out +' '+directory_grid)
        else: 

            print('Already regridded: ' + str(file))

# move grid to regrid folder
os.system('mv '+gridname+".txt "+directory_grid)