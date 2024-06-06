#!/usr/bin/env python3

import rpy2.robjects as robjects
import numpy as np
import os

robjects.r['load']('/mnt/HDS_MEDYCLIM/MEDYCLIM/PREDANTAR/climate_data/obs/oracle/CLIMATE_GRID/coord.lst.Rdata')
matrix1=robjects.r['coord.lst']
coord=np.array(matrix1)
lat0=coord[2,:]
lon0=coord[1,:]

str1="/mnt/HDS_MEDYCLIM/MEDYCLIM/PREDANTAR/climate_data/obs/oracle/CLIMATE_GRID/climate_atlas_"
str2="_CLIMATE_GRID_1950_2022_bis.Rdata"
temp_list=[]
temp_list.append("HUMIDITY_RELATIVE")
temp_list.append("PRECIP_1H_MAX")
temp_list.append("PRECIP_DURATION")
temp_list.append("PRECIP_QUANTITY")
temp_list.append("PRESSURE")
temp_list.append("SHORT_WAVE_FROM_SKY")
temp_list.append("SUN_DURATION")
temp_list.append("SUN_INT")
temp_list.append("TEMP_AVG")
temp_list.append("TEMP_MAX")
temp_list.append("TEMP_MIN")
temp_list.append("WIND_PEAK_SPEED")
temp_list.append("WIND_SPEED")
temp_list.append("SUN_INT_HORIZ")

for variables in temp_list:
  file_name = str1+variables+str2
  robjects.r['load'](file_name)
  matrix3=robjects.r['time.arr']
  pr=np.array(matrix3)
  from netCDF4 import Dataset
  str3 = "/mnt/HDS_FORESTFLOW/FORESTFLOW/CLIMATE_GRID/Daily_"
  str4 = "_1950_2022_Obs_ClimateGrid.nc"
  nc_file_name = str3+variables+str4
  ncfile = Dataset(nc_file_name,mode='w',format='NETCDF4_CLASSIC')
  lat_dim = ncfile.createDimension('lat', 1360)
  lon_dim = ncfile.createDimension('lon', 1360)
  year_dim = ncfile.createDimension('year', 72) 
  month_dim = ncfile.createDimension('month', 12)
  day_dim = ncfile.createDimension('day', 31)    
  lat = ncfile.createVariable('lat', np.float32, ('lat',))
  lat.units = 'Degree North'
  lat.long_name = 'latitude'
  lon = ncfile.createVariable('lon', np.float32, ('lon',))
  lon.units = 'Degree East'
  lon.long_name = 'longitude'
  year = ncfile.createVariable('year', np.float32, ('year',))
  year.units = 'years since 1950'
  year.long_name = 'year'
  month = ncfile.createVariable('month', np.float32, ('month',))
  month.units = 'month'
  month.long_name = 'month'
  day = ncfile.createVariable('day', np.float32, ('day',))
  day.units = 'day'
  day.long_name = 'day'
  precip = ncfile.createVariable(variables,np.float64,('lat','year','month','day'))
  precip.units = ''
  precip.standard_name = variables
  precip[:,:,:,:] = pr
  nmonth = len(month_dim);
  nyear = len(year_dim);
  nday = len(day_dim);
  month[:]=np.arange(nmonth)
  day[:]=np.arange(nday)
  year[:]=np.arange(nyear)
  lat[:]=lat0
  lon[:]=lon0
  ncfile.close()

