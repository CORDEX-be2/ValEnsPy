# Script to download ERA5 or ERA5-Land
"""
Python script to download ERA5 through the CDS API

Necessary dependency cdsapi, see https://cds.climate.copernicus.eu/api-how-to

Adapted for Hortense, direct download to the shared s4c project_input folder of ERA5. Keep the downloads clean!

Resolution of ERA5: 0.25° x 0.25°
Horizontal resolution ERA5-Land: 0.1° x 0.1°

credits: based on the era_single_level by Bert Van Schaeybroeck, edited by Inne Vanderkelen, May 2024

"""


import cdsapi

c = cdsapi.Client()
import os

## USER SETTINGS
dataset = "era5-land"  #'era5' #era5-land' #'cams'
region = "belgium"  # 'europe', 'globe', 'world'; coordinates predefined below.
init_yr = 1940
end_yr = 2024
time_freq = "hourly"  # hourly is base resolution


# set directory
base_dir = "/dodrio/scratch/projects/2022_200/project_input/External/observations/"

var_lst = ["2m_temperature"]


"""
ERA5: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
ERA5-Land: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview
overview of downloadable variables
VARS = [
	'2m_dewpoint_temperature',
	'2m_temperature',
	'skin_temperature',
	'total_precipitation',
	'10m_u_component_of_wind',
	'10m_v_component_of_wind'
	'100m_u_component_of_wind',
	'100m_v_component_of_wind',
	'total_sky_direct_solar_radiation_at_surface',
	'surface_solar_radiation_downward_clear_sky',
	'surface_solar_radiation_downwards',
	'mean_sea_level_pressure',
	'surface_pressure',
	'Boundary_layer_height',
	'skin_temperature',
	'total_column_water_vapour',
	'total_column_water'
	]
snow = [
	'convective_snowfall',
	'convective_snowfall_rate_water_equivalent',
	'large_scale_snowfall',
	'large_scale_snowfall_rate_water_equivalent',
	'snow_albedo',
	'snow_density',
	'snow_depth',
	'snow_evaporation',
	'snowfall',
	'snowmelt',
	'temperature_of_snow_layer',
	'total_column_snow_water'

era5-land = [
	'total_sky_direct_solar_radiation_at_surface'
	'surface_net_solar_radiation',
	'total_evaporation'
	]
cams = [ 'particulate_matter_10um',
	'particulate_matter_1um',
	'particulate_matter_2.5um',
	'total_column_ozone',
	'total_column_nitrogen_dioxide',
	'total_column_nitrogen_monoxide'
	]
"""

times_to_use = [
    "00:00",
    "01:00",
    "02:00",
    "03:00",
    "04:00",
    "05:00",
    "06:00",
    "07:00",
    "08:00",
    "09:00",
    "10:00",
    "11:00",
    "12:00",
    "13:00",
    "14:00",
    "15:00",
    "16:00",
    "17:00",
    "18:00",
    "19:00",
    "20:00",
    "21:00",
    "22:00",
    "23:00",
]

# dataset specific settings (era5, era5-land, cams)
if dataset == "era5":
    grid_to_use = "0.25/0.25"
    dataset_to_use = "reanalysis-era5-single-levels"
    if (init_yr or end_yr) < 1940:
        print("year selected too low, start dataset at 1940")
    if (init_yr or end_yr) > 2024:
        print("year selected in the future, dataset ranges up to 2024")
elif dataset == "era5-land":
    grid_to_use = "0.10/0.10"
    dataset_to_use = "reanalysis-era5-land"

elif dataset == "cams":
    grid_to_use = "0.10/0.10"
    dataset_to_use = "cams-global-reanalysis-eac4"
    times_to_use = (
        ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"],
    )
else:
    print("dataset unknown:", dataset)
    exit()

yr_lst = [x for x in map(str, range(int(init_yr), int(end_yr) + 1))]


# other settings
expver_to_use = "1"
product_type_to_use = "reanalysis"
format_to_use = "netcdf"

# format: lat_min/lon_min/lat_max/lon_max
if region == "europe":
    area_to_use = "Europe"
elif region == "belgium":
    area_to_use = "49./2./52./7."
    grid_to_use = "0.125/0.125"
elif region == "world":
    area_to_use = "-56.5/-179.5/85.5/179.5"
    grid_to_use = "1/1"
elif region == "globe":
    area_to_use = "90/-180/-90/180"
else:
    print("area unknown")
    exit()

# prior to downloading ERA5 data, copy this script into folder.

# download ERA5 data
for var_to_use in var_lst:
    dir_out = (
        base_dir + dataset + "/" + region + "/" + time_freq + "/" + var_to_use + "/"
    )
    print(dir_out)
    os.makedirs(dir_out, exist_ok=True)

    # do actual download
    for yr_to_use in yr_lst:
        print(yr_to_use)
        target_filename = (
            dir_out
            + dataset
            + "-"
            + time_freq
            + "-"
            + region
            + "-"
            + var_to_use
            + "-"
            + yr_to_use
            + ".nc"
        )
        if not os.path.isfile(target_filename):
            print(target_filename)
            c.retrieve(
                dataset_to_use,
                {
                    "product_type": product_type_to_use,
                    "variable": var_to_use,
                    "expver": expver_to_use,
                    "area": area_to_use,
                    "grid": grid_to_use,
                    "year": yr_to_use,
                    "month": [
                        "01",
                        "02",
                        "03",
                        "04",
                        "05",
                        "06",
                        "07",
                        "08",
                        "09",
                        "10",
                        "11",
                        "12",
                    ],
                    "day": [
                        "01",
                        "02",
                        "03",
                        "04",
                        "05",
                        "06",
                        "07",
                        "08",
                        "09",
                        "10",
                        "11",
                        "12",
                        "13",
                        "14",
                        "15",
                        "16",
                        "17",
                        "18",
                        "19",
                        "20",
                        "21",
                        "22",
                        "23",
                        "24",
                        "25",
                        "26",
                        "27",
                        "28",
                        "29",
                        "30",
                        "31",
                    ],
                    "time": times_to_use,
                    "format": format_to_use,
                },
                target_filename,
            )
