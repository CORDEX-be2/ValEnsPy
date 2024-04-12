#Create functionality which checks if the an xarray dataset is a CF compliant dataset
import xarray as xr
from yaml import safe_load
from typing import Union, List

#Get the current file directory and load the variables.yml file
from pathlib import Path
files = Path(__file__).resolve().parent

with open(files / 'data' / 'CORDEX_variables.yml') as file:
    CORDEX_VARIABLES = safe_load(file)

#ToDo add warning messages for the users
def is_cf_compliant(netCDF: Union[str, Path, xr.Dataset]) -> bool:
    if isinstance(netCDF, str) or isinstance(netCDF, Path):
        if _check_file_extension(netCDF):
            try:
                ds = xr.open_dataset(netCDF)
            except:
                return False
        else:
            return False
    else:
        ds = netCDF

    var_meta_data_ok = all([_check_variable_metadata(ds[var]) for var in ds.data_vars if "_bnds" not in var])
    main_meta_data_ok = _check_main_metadata(ds)
    cordex_vars_data_ok = all([_check_variable_by_name(ds[var]) for var in ds.data_vars])
    
    return var_meta_data_ok and main_meta_data_ok and cordex_vars_data_ok

def _check_variable_by_name(da: xr.DataArray):
    """
    Check if the variable is a CORDEX variable. 
    If it is check if the attributes in the datarray match the corresponding attributes in the CORDEX_VARIABLES dictionary

    Parameters
    ----------
    da : xr.DataArray
        The data array to check
    
    Returns
    -------
    bool
        True if the variable is a CORDEX variable and the attributes match or if the variable is not a CORDEX variable, False otherwise
    """
    if da.name in CORDEX_VARIABLES: 
        return all([CORDEX_VARIABLES[da.name][attr] == da.attrs[attr] for attr in da.attrs if attr in CORDEX_VARIABLES[da.name]])
    else:
        return True

def _check_main_metadata(ds: xr.Dataset):
    """
    Check if the required main metadata attributes exist

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to check
    
    Returns
    -------
    bool
        True if the required main metadata attributes exist, False otherwise
    """

    MAIN_METADATA=['title', 'history']

    return all([attr in ds.attrs for attr in MAIN_METADATA])

def _check_variable_metadata(da: xr.DataArray):
    """
    Check if the CF required variable attributes exist

    Parameters
    ----------
    da : xr.DataArray
        The data array to check
    
    Returns
    -------
    bool
        True if the CF required variable attributes exist, False otherwise
    """

    VARIABLE_METADATA = ['units', 'standard_name', 'long_name']
    return all([attr in da.attrs for attr in VARIABLE_METADATA]) 

def _check_file_extension(file: Union[str, Path]) -> bool:
    """Check if the file extension is netcdf (.nc)"""
    return file.endswith('.nc')

if __name__ == '__main__':
    print(is_cf_compliant("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/tests/data/tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc"))
    EOBS_data_dir = Path("/dodrio/scratch/projects/2022_200/project_input/External/observations/EOBS/0.1deg/")

    EOBS_obs_files = list(EOBS_data_dir.glob("*tg*mean*.nc")) #Select all the netCDF files in the directory

    EOBS_ds = xr.open_mfdataset(EOBS_obs_files, combine='by_coords', chunks='auto')
    print(is_cf_compliant(EOBS_ds))

    EOBS_ds.attrs['title'] = 'EOBS dataset'
    print(is_cf_compliant(EOBS_ds))

    EOBS_ds = EOBS_ds.rename({'tg': 'tas'})
    print(is_cf_compliant(EOBS_ds))