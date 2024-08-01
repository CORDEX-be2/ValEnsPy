# Helper functions for the regridding

import numpy as np
import xarray as xr
import valenspy as vp
from valenspy._utilities import  load_yml
from valenspy.preprocessing_tasks.regrid import remap_cdo
from pathlib import Path
import os
 

 # general settings
machine = "hortense"
src_path = Path("../src/valenspy")
model = 'CCLM'
postproc_base_dir = "/dodrio/scratch/projects/2022_200/RCS/CORDEXBE2/postprocessing/"
# get lookup file for model
mod_LOOKUP = load_yml(model+"_lookup")
    
# define bounds 
bounds = {'europe':
                    {'lat_bounds': [70,35], 
                    'lon_bounds': [-15,40]}, 
        'belgium': 
                    {'lat_bounds': [52,49], 
                    'lon_bounds': [2,7]}}

def remap_cdo_intermediatefiles(gridfile: str, ref_dataset: str, files_to_regrid: list, remap_method: str):
    """
    Remaps a list of NetCDF files to a new grid using CDO (Climate Data Operators).

    This function performs regridding of the given NetCDF files to the grid defined 
    by `gridfile` using the specified remap method. It saves the regridded files in 
    a "regridded" subdirectory within the directory where the original files are located 
    and returns the combined dataset.

    Parameters:
    -----------
    gridfile : str
        Path to the NetCDF file that defines the target grid for regridding.

    ref_dataset : str
        Reference dataset name to be used in the naming of the regridded files.
    
    files_to_regrid : list
        List of paths to the NetCDF files that need to be regridded.
    
    remap_method : str
        Method to be used for regridding. This should be a valid remapping method 
        supported by CDO (e.g., "remapbil", "remapcon").

    Returns:
    --------
    xarray.Dataset
        The combined xarray Dataset containing all the regridded files.

    Notes:
    ------
    - This function uses the CDO command line tool to perform the regridding. 
      Ensure that CDO is installed and accessible from the command line.
    - If using the "remapcon" method, ensure that the grid cell bounds are provided 
      in the input files. Additional processing might be required to add these bounds.
    - The function renames longitude and latitude dimensions to "lon" and "lat" 
      respectively, if they are named differently in the regridded files.
    """
    # within the directory where the files to be regridded are stored, create a "regridded" folder to store the regridded files
    regrid_dir = os.path.dirname(files_to_regrid[0]) + "/regridded/"

    if not os.path.exists(regrid_dir):
        os.makedirs(regrid_dir)

    for file in files_to_regrid:
        regridded_file = os.path.join(regrid_dir, os.path.basename(file)[:-3] + '_' + ref_dataset + '_grid.nc')
        cdo_string = f"cdo {remap_method},{gridfile} {file} {regridded_file}"
        os.system(cdo_string)

    files_regridded_all = os.path.join(regrid_dir, '*' + '_' + ref_dataset + '_grid.nc')
    ds = xr.open_mfdataset(files_regridded_all, combine="by_coords", chunks="auto")

    # Rename dimensions if necessary
    ds = ds.rename({"longitude": "lon", "latitude": "lat"}) if "longitude" in ds.coords or "latitude" in ds.coords else ds

    return ds




import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def calc_and_plot_timmean_bias(region: str, bounds:dict, da_mod, da_obs, experiment, ref_dataset): 

    lat_bounds = bounds[region]['lat_bounds']
    lon_bounds = bounds[region]['lon_bounds']

    da_mod = da_mod.sel(lon=slice(lon_bounds[0],lon_bounds[1]),lat=slice(lat_bounds[0], lat_bounds[1])).mean('time', keep_attrs=True)
    da_obs = da_obs.sel(lon=slice(lon_bounds[0],lon_bounds[1]),lat=slice(lat_bounds[0], lat_bounds[1])).mean('time', keep_attrs=True)
    bias = da_mod - da_obs


    fig, axes = plt.subplots(1, 3, figsize=(14, 3), subplot_kw={'projection': ccrs.PlateCarree()})
    axes = axes.flatten()

    # find plotting min and max
    vmin = float(min(da_mod.min().values, da_obs.min().values))
    vmax = float(max(da_mod.max().values, da_obs.max().values))

    bias_bound = float(max(abs(bias.min().values), abs(bias.max().values)))


    cbar_label = f"{da_obs.attrs['long_name']} ({da_obs.attrs['units']})"
    cbar_kwargs = {'label': cbar_label}


    # Define a function to add borders, coastlines to the axes
    def add_features(ax):
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.COASTLINE, linewidth = 0.5, color = 'k')
        ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())

    # obs
    ax = axes[0]
    da_obs.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
    ax.set_title(da_obs.attrs['long_name'], loc='left')
    ax.set_title('')
    ax.set_title(ref_dataset, loc='right')
    add_features(ax)

    # mod
    ax = axes[1]
    da_mod.plot(ax=ax, vmin=vmin, vmax=vmax, cbar_kwargs=cbar_kwargs)
    ax.set_title('')
    ax.set_title(experiment, loc='right')
    add_features(ax)

    # bias
    ax = axes[2]
    bias.plot(ax=ax, cmap = 'coolwarm', vmax = bias_bound, vmin = - bias_bound, cbar_kwargs=cbar_kwargs)
    ax.set_title('')
    ax.set_title(f"{experiment} - {ref_dataset}", loc='right')
    add_features(ax)

    fig.tight_layout()
    plt.show()

    return fig
    #fig.savefig(f"./plots/{region}_{experiment}_{ref_dataset}_timmean_bias.png")


def geo_to_rot(coord, ds):
    ###
    # Converts a geographic (lon,lat) point to an (rlon,rlat) point
    # => coord: geographic lon - lat couple (degrees) => list or tuple
    # => ds: xarray dataset of a CCLM file
    ###
    # Read in COSMO-CLM data

    rp_lat = float(ds.rotated_pole.grid_north_pole_latitude)
    rp_lon = float(ds.rotated_pole.grid_north_pole_longitude)
    # Convert 
    co = np.deg2rad(coord)
    rp_lat = np.deg2rad(rp_lat); rp_lon = np.deg2rad(rp_lon)
    p_rlat = np.arcsin(np.sin(co[1])*np.sin(rp_lat) + np.cos(co[1])*np.cos(rp_lat)*np.cos(co[0]-rp_lon)) 
    p_rlon = np.arctan((np.cos(co[1])*np.sin(co[0]-rp_lon)) / (np.cos(co[1])*np.sin(rp_lat)*np.cos(co[0]-rp_lon) - np.sin(co[1])*np.cos(rp_lat))) 
    p_rlat = np.rad2deg(p_rlat); p_rlon = np.rad2deg(p_rlon)
    # Return rlon-rlat couple
    return [p_rlon,p_rlat]
    



def load_calc_plot_bias_map(variable: str, ref_dataset: str, experiments: list, months_to_analyse: list, region='europe', unit_conversion = False): 

    # ------------------------------
    # 1. Load reference data

    # start up input manager
    manager = vp.InputManager(machine=machine)

    # use input manager to load data, defined on settings above
    ds_obs = manager.load_data(ref_dataset,variable, period=[1995,1995],freq="hourly",region=region)
    ds_obs = ds_obs.resample(time='1D').mean()    

    # quick and dirty fix to account for correct units for cumulative variables in ERA5 and ERA5 land - to be properly solved in the inputmanager
    if unit_conversion == True: 
        print('did unit conversion')
        ds_obs[variable] = ds_obs[variable]/(86400)
    elif unit_conversion=='pr': 
        print('did pr unit conversion')
        ds_obs[variable] = ds_obs[variable]/(86400)*1000
    # retrieve ERA5 gridfile - for regridding 
    gridfile = manager._get_file_paths(ref_dataset,variable, period=[1995,1995],freq="hourly",region=region)[0]


    # ------------------------------
    # 2. Load and regrid model data

    for experiment in experiments: 

        print(experiment)

        mod_LOOKUP = load_yml(model+"_lookup")
        # get CCLM variable corresponding to the requested variable using its look-up table
        mod_var = mod_LOOKUP[variable]['mod_name']

        # define the path
        directory = Path(postproc_base_dir + experiment +'/'+mod_var + '/')

        # define the CCLM files for the corresponding variable
        mod_files = list(directory.glob(mod_var+"_daymean.nc")) # Select all the netCDF files in the directory

        if not mod_files:  # empty list
                print(f"{variable} not available for {experiment}")
        else: 
            
            ds_mod = xr.open_mfdataset(mod_files, combine="by_coords", chunks="auto")

            # regrid
            ds_mod = remap_cdo_intermediatefiles(gridfile, ref_dataset,  mod_files, remap_method = "remapbil")
            # ds_mod = remap_cdo(gridfile, ds_mod, remap_method = "con")
            # ------------------------------
            # 3. Calculate diagnostic and do plotting
            
            # select variable and corresponding period
            da_mod = ds_mod[mod_var].sel(time=ds_mod['time'].dt.month.isin(months_to_analyse))
            da_obs = ds_obs[variable].sel(time=ds_obs['time'].dt.month.isin(months_to_analyse))

            calc_and_plot_timmean_bias(region, bounds, da_mod, da_obs, experiment, ref_dataset)


# ------------------------------
# 1. Load reference data


def plot_point_timeseries(variable: str, ref_dataset: str, experiments: list, point_coord: tuple, point_id:str,  months_to_analyse: list): 
        
    # start up input manager
    manager = vp.InputManager(machine=machine)

    # use input manager to load data, defined on settings above
    ds_obs = manager.load_data(ref_dataset,variable, period=[1995,1995],freq="daily",region=region, path_identifiers = ["-daily-"])

    # select point 
    ds_obs_point = ds_obs.sel(lon=point_coord[0],lat=point_coord[1], method='nearest')
    da_obs_point = ds_obs_point[variable].sel(time=ds_obs['time'].dt.month.isin(months_to_analyse))


    # generate path of CCLM output

    # dictorionary to save data arrays of experiments
    d_da_mod_point = {}

    for experiment in experiments: 
        # get CCLM variable corresponding to the requested variable using its look-up table
        mod_var = mod_LOOKUP[variable]['mod_name']

        # define the path
        directory = Path(postproc_base_dir + experiment +'/'+mod_var + '/')

        # open the CCLM file for the corresponding variable
        mod_files = list(directory.glob(mod_var+"_daymean.nc")) # Select all the netCDF files in the directory

        if not mod_files:  # empty list - move to next element in loop
            print(f"{variable} not available for {experiment}")
            continue
            
        ds_mod = xr.open_mfdataset(mod_files, combine="by_coords", chunks="auto")

        coord_points_rotated = geo_to_rot(point_coord, ds_mod)
        ds_mod_point = ds_mod.sel(rlon = coord_points_rotated[0], rlat = coord_points_rotated[1], method='nearest')

        da_mod_point = ds_mod_point[mod_var].sel(time=ds_mod['time'].dt.month.isin(months_to_analyse))

        d_da_mod_point[experiment] = da_mod_point

    # do plotting
    fig, ax = plt.subplots(figsize = (7,3))

    da_obs_point.plot(ax=ax, label = ref_dataset, color='k')

    for experiment in d_da_mod_point: 
        d_da_mod_point[experiment].plot(ax=ax, label = experiment, alpha=0.5)

    ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_title(da_obs_point.attrs['long_name'], loc='left')
    ax.set_title(' ', loc='center')
    ax.set_title(f"{point_id} ({point_coord[1]}°N,  {point_coord[0]}°E)", loc='right'); 

    fig.savefig(f"./plots/timeseries_{variable}_{point_coord[1]}N_{point_coord[0]}E).png")