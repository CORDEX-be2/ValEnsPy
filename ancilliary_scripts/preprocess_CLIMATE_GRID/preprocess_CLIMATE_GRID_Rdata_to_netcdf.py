#!/usr/bin/env python3

# Script to load the CLIMATE DATA .Rdata files and grid them and convert them to netcdf files per variable
# uses the .RData files retrieved from the RMI Oracle databse through the climate_grid_daily.txt script from Bert Van Schaeybroeck

# I. Vanderkelen, June 2024

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import xarray as xr
import rioxarray
import pandas as pd
from datetime import date
import numpy as np
import cartopy.crs as ccrs


# activate r to pandas convertor
pandas2ri.activate()


# user settings
variables = ["EVAPOTRANS_REF", "SUN_INT", "SUN_DURATION", "PRECIP_DURATION", "WIND_PEAK_SPEED", "PRECIP_1H_MAX", "EVAPOTRANS_REF", "TEMP_MAX","HUMIDITY_RELATIVE","TEMP_MIN", "TEMP_AVG", "WIND_SPEED", "PRESSURE", "SHORT_WAVE_FROM_SKY", "SUN_INT_HORIZ", "PRECIP_QUANTITY"]


# INFO ON GRID
# load climategrid meta data on variables and units
meta = pd.read_csv('CLIMATE_GRID_meta.csv', delimiter=";")

# load information on projection used 
proj_string = "+proj=lcc +lat_1=49.83333388888889 +lat_2=51.16666722222222 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +units=m +no_defs"
# load the pixel lat and lon variable and use this to transpose to own defined grid
df_coords_points = pd.read_csv('grid_5kmx5km.csv', header=1, delimiter=' ') # these are the lat lons and lambert coordinates for all pixels in CLIMATE_DATA


# load the full grid, creased based on the proj_string and following bounding points:
#NE_lon, NE_lat = 9.53269211610237, 53.4367017362904
#SW_lon, SW_lat = 0.163155782953472, 47.515819098539
df_full_grid = pd.read_csv('lambert_coordinates_full_climate_grid.csv') # made using R script Michel Journee
lambert_x_grid_raw = df_full_grid['x1'].unique()
lambert_y_grid_raw = df_full_grid['x2'].unique()

# cut this grid to boundingbox including gridcells from CLIMATE_GRID
lambert_x_grid_cutlow = lambert_x_grid_raw[lambert_x_grid_raw>=df_coords_points['LAMBERT_X'].min()]
lambert_x_grid = lambert_x_grid_cutlow[lambert_x_grid_cutlow<=  df_coords_points['LAMBERT_X'].max()]
lambert_y_grid_cutlow = lambert_y_grid_raw[lambert_y_grid_raw>=df_coords_points['LAMBERT_Y'].min()]
lambert_y_grid = lambert_y_grid_cutlow[lambert_y_grid_cutlow<=  df_coords_points['LAMBERT_Y'].max()]


# Find the nearest index in the lons and lats grids and add this to coordinates dataframe
def find_nearest(array, value):
	idx = (np.abs(array - value)).argmin()
	return idx

df_coords_points['LAMBERT_X_INDEX'] = df_coords_points['LAMBERT_X'].apply(lambda x: find_nearest(lambert_x_grid, x))
df_coords_points['LAMBERT_Y_INDEX'] = df_coords_points['LAMBERT_Y'].apply(lambda x: find_nearest(lambert_y_grid, x))


data_dir = '/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc31332_inne/CLIMATE_GRID/' # datadir on hortense
data_dir = '/mnt/HDS_CLIMATE/CLIMATE/CLIMATE_GRID/' # datadir on kili

for variable in variables: 
	print('Converting '+ variable)

	filename = 'climate_atlas_'+str(variable)+'_CLIMATE_GRID_1950_2023.Rdata'

	# load the robject
	robjects.r['load'](data_dir+filename)


	# load the time data and convert to dates
	time = robjects.r['time.lst']
	dates = pd.to_datetime(time)

	# load the vector data
	data = robjects.r['grid.vec']

	# create empty array to fill with gridded data - also for lat and lon
	grid_data = np.full(( len(dates),  len(lambert_y_grid),len(lambert_x_grid)  ), np.nan)

	lat_2d = np.full((len(lambert_y_grid),len(lambert_x_grid)  ), np.nan)
	lon_2d = np.full((len(lambert_y_grid),len(lambert_x_grid)  ), np.nan)

	# Fill the grid data array
	for index, row in df_coords_points.iterrows():
		lambert_x_idx = int(row['LAMBERT_Y_INDEX'])
		lambert_y_idx = int(row['LAMBERT_X_INDEX'])

		pixel_id = int(row['PIXEL_ID'])
		
		grid_data[:, lambert_x_idx, lambert_y_idx ] = data[int(pixel_id) - 1, :]

		lat_2d[ lambert_x_idx, lambert_y_idx ] = df_coords_points[df_coords_points['PIXEL_ID'] == pixel_id]["LAT"].values[0]
		lon_2d[ lambert_x_idx, lambert_y_idx ] = df_coords_points[df_coords_points['PIXEL_ID'] == pixel_id]["LON"].values[0]


	# get metadata from meta dataframe
	unit = meta.loc[meta['variable'] == variable, 'unit'].values[0] 
	long_name = meta.loc[meta['variable'] == variable, 'long_name'].values[0] 
	description = meta.loc[meta['variable'] == variable, 'description'].values[0] 

	# create data array
	da = xr.DataArray(
		data=grid_data,
		dims=["time", "y", "x"],
		coords=dict(
			y=lambert_y_grid,
			x=lambert_x_grid,
			time=dates,
		),
		attrs=dict(
			long_name=long_name,
			description = description,
			units=unit,
		),
	)


	da['x'].attrs = {'units':"E[east]: Easting (meters)", 'long_name': " x coordinate Lambert Conic Conformal (2SP)"}
	da['y'].attrs = {'units':"N[north]: Northing (meters)", 'long_name': "y coordinate Lambert Conic Conformal (2SP)"}

	# convert to dataset and give dataset attributes
	ds = da.to_dataset(name=variable)

	# add 2d lat and lon in lambert coordinates as well
	ds["lat"] = xr.DataArray(
		data=lat_2d,
		dims=["y", "x"],
		coords=dict(
			y=lambert_y_grid,
			x=lambert_x_grid,
		),
		attrs=dict(
			long_name="latitude",
			description = "WGS84 latitude, from values of CLIMATE_GRID, provided per grid point",
			units="degrees_north",
		),
	)
	
	# interpolatie to also have lat values outside of Belgium
	ds['lat'] = ds['lat'].interpolate_na(dim='x', method='linear', fill_value="extrapolate")

	ds["lon"] = xr.DataArray(
		data=lon_2d,
		dims=["y", "x"],
		coords=dict(
			y=lambert_y_grid,
			x=lambert_x_grid,
		),
		attrs=dict(
			long_name="longitude",
			description = "WGS84 longitude, from values of CLIMATE_GRID, provided per grid point",
			units="degrees_east",
		),
	)
	d_attrs = {"creation_date": date.today().strftime("%d-%m-%Y"),
	"creators": "Ghilain N., Van Schaeybroeck B., Vanderkelen I.", 
	"contact": "inne.vanderkelen@meteo.be",
	"version": "1.1", "affiliation": "Royal Meteorological Institute of Belgium", 
	"projection":proj_string}
	ds.attrs = d_attrs

	# interpolatie to also have lon values outside of Belgium
	ds['lon'] = ds['lon'].interpolate_na(dim='x', method='linear', fill_value="extrapolate")

	# also pass crs using rioxarray - passing the spatial_ref
	ds.rio.write_crs(ccrs.Projection(proj_string), inplace=True)

	# export to netcdf
	filename_out = str(variable)+'_CLIMATE_GRID_'+str(dates.year.min())+'_'+str(dates.year.max())+'_daily.nc'
	ds.to_netcdf(data_dir + filename_out, encoding={'time':  {'dtype': 'int32'} })

