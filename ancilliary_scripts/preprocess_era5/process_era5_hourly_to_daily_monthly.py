# Script to merge hourly ERA5 or ERA5-Land into daily or monthly

"""
credits: based on the step1_convert_era5_to_daily_monthly.txt R script by Bert Van Schaeybroeck, edited by Inne Vanderkelen, August 2024
"""
import os
import subprocess
import shutil

import pandas as pd


# Directories setup
dir_base = "/dodrio/scratch/projects/2022_200/project_input/External/observations/"

# Variables initialization
flag_calc_min_max = False # if set to true, also the min and max are calculated, in addition to the mean. 
dataset = "era5"  # "era5-land"
domain = "europe"  # "europe"
init_yr = 1940
end_yr = 2023

var_name_lst =  ["2m_temperature"]

# all possible variables
"""
era5 variables
 variables = [
    "100m_u_component_of_wind", "Boundary_layer_height", "soil_temperature_level_1", "total_cloud_cover", 
    "total_column_water_vapour", "volumetric_soil_water_layer_2", "100m_v_component_of_wind", "high_cloud_cover", 
    "soil_temperature_level_2", "total_column_cloud_ice_water", "total_precipitation", "volumetric_soil_water_layer_3", 
    "10m_u_component_of_wind", "low_cloud_cover", "soil_temperature_level_3", "total_column_cloud_liquid_water", 
    "total_sky_direct_solar_radiation_at_surface", "volumetric_soil_water_layer_4", "10m_v_component_of_wind", 
    "mean_sea_level_pressure", "soil_temperature_level_4", "total_column_rain_water", "u_component_of_wind-500hPa", 
    "wa500", "2m_dewpoint_temperature", "medium_cloud_cover", "soil_water_level_1", "total_column_snow_water", 
    "v_component_of_wind-500hPa", "was", "2m_temperature", "skin_temperature", "surface_solar_radiation_downwards", 
    "total_column_water", "volumetric_soil_water_layer_1"
]
"""
"""
era_5-land variables
variables = [
    "10m_u_component_of_wind", "potential_evaporation", "soil_temperature_level_1", "Surface_net_thermal_radiation", 
    "total_precipitation", "10m_v_component_of_wind", "runoff", "soil_temperature_level_2", "surface_pressure", 
    "volumetric_soil_water_layer_1", "2m_dewpoint_temperature", "skin_temperature", "soil_temperature_level_3", 
    "surface_runoff", "volumetric_soil_water_layer_2", "2m_temperature", "snow_depth", "soil_temperature_level_4", 
    "surface_sensible_heat_flux", "volumetric_soil_water_layer_3", "evaporation", "snow_depth_water_equivalent", 
    "surface_latent_heat_flux", "surface_solar_radiation_downwards", "volumetric_soil_water_layer_4", 
    "mean_sea_level_pressure", "snowfall", "surface_net_solar_radiation", "surface_thermal_radiation_downwards"
]
"""



if flag_calc_min_max:
    statistic = ["mean", "max", "min"]
    statistic_text = {"mean": "", "max": "_max", "min": "_min"}
else:
    statistic = ["mean"]
    statistic_text = {"mean": ""}

amt_fun = len(statistic)

amt_var = len(var_name_lst)
yr_lst = range(init_yr, end_yr + 1)
amt_yr = len(yr_lst)

time_freq_hr = "hourly"  # Or whichever time frequency is applicable
time_freq_day = "daily"
time_freq_mon = "monthly"


for var_to_use in var_name_lst:
    print(f"dataset: {dataset}, var: {var_to_use}")
    dir_base_var_hr = os.path.join(dir_base, dataset, domain, var_to_use, time_freq_hr) 
    dir_base_var_day = os.path.join(dir_base, dataset, domain, var_to_use, time_freq_day) 
    os.makedirs(dir_base_var_day, exist_ok=True)
    dir_base_var_mon = os.path.join(dir_base, dataset, domain, var_to_use, time_freq_mon) 
    os.makedirs(dir_base_var_mon, exist_ok=True)
    if dataset == "era5-land": 
        df_vars = pd.read_csv('era5-land_vars.csv')
        flag_cumul = df_vars.loc[df_vars['name.long'] == var_to_use,'cumul'] ==1
    else: 
        flag_cumul = False
    for yr_to_use in yr_lst:
        print(f"var: {var_to_use}, year: {yr_to_use}")
        
        file_base_hr = f"{dataset}-{time_freq_hr}-{domain}-{var_to_use}-{yr_to_use}.nc"
        file_hr = os.path.join(dir_base_var_hr, file_base_hr)
        
        if flag_cumul:
            file_tmp1 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp1")
            file_tmp2 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp2")
            file_tmp3 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp3")
            subprocess.run(f"cdo -b F64 -shifttime,-1hour -selhour,1 {file_hr} {file_tmp1}", shell=True)
            subprocess.run(f"cdo -b F64 -shifttime,-1hour -selhour,0,2/23 -deltat {file_hr} {file_tmp2}", shell=True)
            subprocess.run(f"cdo -mergetime {file_tmp1} {file_tmp2} {file_tmp3}", shell=True)
            os.remove(file_tmp1)
            os.remove(file_tmp2)
        elif var_to_use == "10m_v_component_of_wind":
            file_tmp1 = file_hr.replace("_v_", "_u_")
            file_tmp2 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp2")
            file_tmp3 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp3")
            subprocess.run(f"cdo merge {file_hr} {file_tmp1} {file_tmp2}", shell=True)
            subprocess.run(f"cdo expr,'v10=sqrt(u10*u10+v10*v10)' {file_tmp2} {file_tmp3}", shell=True)
            os.remove(file_tmp2)
        elif var_to_use == "10m_u_component_of_wind":
            file_tmp1 = file_hr.replace("_u_", "_v_")
            file_tmp2 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp2")
            file_tmp3 = os.path.join(scratch_dir_to_use, f"{file_base_hr}.tmp3")
            subprocess.run(f"cdo merge {file_hr} {file_tmp1} {file_tmp2}", shell=True)
            subprocess.run(f"cdo expr,'u10=sqrt(u10*u10+v10*v10)' {file_tmp2} {file_tmp3}", shell=True)
            os.remove(file_tmp2)
        else:
            file_tmp3 = file_hr

        for fun_to_use in statistic:
            file_day = os.path.join(dir_base_var_day, f"{dataset}-{time_freq_day}-{domain}-{var_to_use}{statistic_text[fun_to_use]}-{yr_to_use}.nc")
            file_mon = os.path.join(dir_base_var_mon, f"{dataset}-{time_freq_mon}-{domain}-{var_to_use}{statistic_text[fun_to_use]}-{yr_to_use}.nc")
            
            if os.path.exists(file_day):
                print(f"{file_day}  exists.")
            else:
                subprocess.run(f"cdo -z zip day{fun_to_use} {file_tmp3} {file_day}", shell=True)

            if os.path.exists(file_mon):
                print(f"{file_mon}  exists.")
            else:
                subprocess.run(f"cdo -z zip monmean {file_day} {file_mon}", shell=True)

        if flag_cumul:
            os.remove(file_tmp3)



