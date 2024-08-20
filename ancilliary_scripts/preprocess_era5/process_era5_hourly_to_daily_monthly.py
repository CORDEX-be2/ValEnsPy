# Script to merge hourly ERA5 or ERA5-Land into daily or monthly

"""
credits: based on the step1_convert_era5_to_daily_monthly.txt by Bert Van Schaeybroeck, edited by Inne Vanderkelen, August 2024
"""

import os
import subprocess
import shutil


# Load required module
subprocess.run(["module", "load", "CDO"], shell=True)

# Directories setup
base_dir = "/dodrio/scratch/projects/2022_200/project_input/External/observations/"

# Variables initialization
flag_calc_min_max = False # if set to true, also the min and max are calculated, in addition to the mean. 
dataset = "era5"  # "era5-land"
domain = "europe"  # "europe"
init_yr = 1940
end_yr = 2023

var_name_lst = [
    "2m_temperature"
]

if min_max_msk:
    statistic = ["mean", "max", "min"]
    statistic_text = {"mean": "", "max": "_max", "min": "_min"}
else:
    statistic = ["mean"]
    statistic_lst = {"mean": ""}

amt_fun = len(fun_lst)

amt_var = len(var_name_lst)
yr_lst = range(init_yr, end_yr + 1)
amt_yr = len(yr_lst)

dir_base = base_dir+dataset+'/'  # Set this to the appropriate base directory
time_freq_hr = "hourly"  # Or whichever time frequency is applicable
time_freq_day = "daily"
time_freq_mon = "monthly"

dir_base_hr = os.path.join(dir_base, dataset, domain, time_freq_hr)
dir_base_day = os.path.join(dir_base, dataset, domain, time_freq_day)
dir_base_mon = os.path.join(dir_base, dataset, domain, time_freq_mon)

# Processing loop
for var_to_use in var_name_lst:
    print(f"dataset: {dataset}, var: {var_to_use}")
    dir_base_var_hr = os.path.join(dir_base_hr, var_to_use)
    dir_base_var_day = os.path.join(dir_base_day, var_to_use)
    os.makedirs(dir_base_var_day, exist_ok=True)
    dir_base_var_mon = os.path.join(dir_base_mon, var_to_use)
    os.makedirs(dir_base_var_mon, exist_ok=True)
    
    use_cumul = False  # This needs to be set based on the specific dataset
    for yr_to_use in yr_lst:
        print(f"var: {var_to_use}, year: {yr_to_use}")
        
        file_base_hr = f"{dataset}-{time_freq_hr}-{domain}-{var_to_use}-{yr_to_use}.nc"
        file_hr = os.path.join(dir_base_var_hr, file_base_hr)
        
        if use_cumul:
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

        for fun_to_use in fun_lst:
            file_day = os.path.join(dir_base_var_day, f"{dataset}-{time_freq_day}{fun_str_lst[fun_to_use]}-{domain}-{var_to_use}-{yr_to_use}.nc")
            file_mon = os.path.join(dir_base_var_mon, f"{dataset}-{time_freq_mon}{fun_str_lst[fun_to_use]}-{domain}-{var_to_use}-{yr_to_use}.nc")

            subprocess.run(f"cdo -z zip day{fun_to_use} {file_tmp3} {file_day}", shell=True)
            subprocess.run(f"cdo -z zip monmean {file_day} {file_mon}", shell=True)

        if use_cumul:
            os.remove(file_tmp3)


