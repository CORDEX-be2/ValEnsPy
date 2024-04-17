#Create functionality which checks if the an xarray dataset is a CF compliant dataset
import xarray as xr
from yaml import safe_load
from typing import Union, List

#Get the current file directory and load the variables.yml file
from pathlib import Path
files = Path(__file__).resolve().parent

with open(files / 'data' / 'CORDEX_variables.yml') as file:
    CORDEX_VARIABLES = safe_load(file)

MAIN_METADATA=['Conventions', 'history']
VARIABLE_METADATA = ['units', 'standard_name', 'long_name']

#Currently only a True or False: True if the file is CF compliant, False otherwise
#Should be extended to return a list of errors if the file is not CF compliant. Maybe each check should through warning separately?
def is_cf_compliant(netCDF: Union[str, Path, xr.Dataset]) -> bool:
    """
    Check if a file is a CF compliant netCDF file. The following checks are performed:
    - Check if the file is a netCDF file (or an xarray dataset)
    - Check if the main metadata attributes exist (title, history)
    - Check if the variable metadata attributes exist (units, standard_name, long_name)
    - For each variable check if it is present in the predefined variables, if so check if the attributes match.

    Parameters
    ----------
    netCDF : Union[str, Path, xr.Dataset]
        The netCDF file to check or the xarray dataset to check
    
    Returns
    -------
    bool
        True if the file is CF compliant, False otherwise

    Examples
    --------
    .. code-block:: python

        >>> import valenspy as vp
        >>>
        >>> vp.cf_checks.is_cf_compliant(vp.demo_data_CF)
        True
    
    """
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

    return all([attr in da.attrs for attr in VARIABLE_METADATA]) 

def _check_file_extension(file: Union[str, Path]) -> bool:
    """Check if the file extension is netcdf (.nc)"""
    if isinstance(file, str):
        return file.endswith('.nc')
    else:
        return file.suffix == '.nc'