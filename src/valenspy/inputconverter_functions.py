from pathlib import Path
import xarray as xr
from yaml import safe_load
from valenspy.cf_checks import is_cf_compliant, cf_status

# get path of source code (current path)
src_path = Path(__file__).resolve().parent 

# open CORDEX variable lookup dictionary
with open(src_path / "ancilliary_data" / "CORDEX_variables.yml") as file:
    CORDEX_VARIABLES = safe_load(file)



def _non_convertor(file: Path) -> Path:
    """A dummy function that does not convert the input file as it already is in the correct format."""
    return file



def EOBS_to_CF(file: Path) -> Path:

    """
    Convert the EOBS netCDF file to a netCDF file in CF convention 

    Parameters
    ----------
    file : Path
        The path to the netCDF file of specific variable to convert 
    
    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """

    # open the observation dataset
    ds = xr.open_mfdataset(file, combine='by_coords', chunks='auto')

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator? 
    obsdata_name = "EOBS"

    with open(src_path / "ancilliary_data" / Path(obsdata_name+"_lookup.yml")) as file:
        obs_LOOKUP = safe_load(file)

    
    # make EOBS CF compliant

    for var_obs in ds.data_vars: 

        # Get the CORDEX variable in the observational dataset using the observational lookup table  
        var = next((k for k, v in obs_LOOKUP.items() if v.get('obs_name') == var_obs), None)

        if var: #Dont processes variables that are not in the lookup table.

            # update variable name to CORDEX variable name
            ds = ds.rename_vars({obs_LOOKUP[var]["obs_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded EOBS units for units different to CORDEX
            if obs_LOOKUP[var]['obs_units'] == 'Celcius': 
                ds[var] = _convert_Celcius_to_Kelvin(ds[var]) 
            elif obs_LOOKUP[var]['obs_units'] == 'hPa': 
                ds[var] = _convert_hPa_to_Pa(ds[var]) # hPa to Pa
            elif obs_LOOKUP[var]['obs_units'] == 'mm': # ! note observations remain daily time frequency
                ds[var] = _convert_mm_to_kg_m2_s1(ds[var]) # mm to kg m^-2 s^-1 

            # update unit attribute
            ds[var].attrs["units"] = CORDEX_VARIABLES[var]["units"] # from the CORDEX look-up table 

            # add necessary metadata
            ds[var].attrs["standard_name"]      = CORDEX_VARIABLES[var]["standard_name"]  # from the CORDEX look-up table
            ds[var].attrs["long_name"]          = CORDEX_VARIABLES[var]["long_name"]  # from the CORDEX look-up table
            ds[var].attrs["original_name"]      = obs_LOOKUP[var]["obs_name"]
            ds[var].attrs["original_long_name"] = obs_LOOKUP[var]["obs_long_name"]

            # rename dimensions
            ds[var] = ds[var].rename({"latitude": "lat", "longitude": "lon"})

            # convert the time dimension to a pandas datetime index --  do we want this to happen within the convertor? Or do we leave it up to the user? 
            ds[var]['time'] = pd.to_datetime(ds[var].time)


            # additional attributes -- hard coded for EOBS
            ds[var].attrs["freq"]               = "daily"  # possible values: daily, hourly, monthly, yearly
            ds[var].attrs["spatial_resolution"] = "0.1deg"  
            ds[var].attrs["domain"]             = "europe"

    
    # Soft check for CF compliance 
    cf_status(ds)

    return ds

def ERA5_to_CF(file: Path) -> Path:

    """
    Convert the ERA5 netCDF file to a xarray Dataset in CF convention 

    Parameters
    ----------
    file : Path
        The path to the netCDF file of specific variable to convert 
    
    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """

    # open the observation dataset
    ds = xr.open_mfdataset(file, combine='by_coords', chunks='auto')

    # Soft check for CF compliance 
    cf_status(ds)
    
    return ds




# helper functions for unit conversion - can be moved to more appropriate place 

## imports for helpen functions
import xarray as xr
import pandas as pd
import numpy as np

# Do we want other possible inputs than data arrays? 
def _convert_Celcius_to_Kelvin(da: xr.DataArray): 

    """
    Convert values in xarray DataArray from °C to K

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """
    
    # do conversion
    da = da + 273.15 # Celcius to Kelvin

    # update units attribute --  naming of units defined in ./src/valenspy/ancilliary_data/CORDEX_variables.yml
    da["units"] = "K" 

    return da

def _convert_Kelvin_to_Celcius(da: xr.DataArray): 

    """
    Convert values in xarray DataArray from K to °C

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """
    
    # do conversion
    da = da - 273.15 # Kelvin to Celcius

    # update units attribute
    da["units"] = "°C" 

    return da

def _convert_hPa_to_Pa(da: xr.DataArray): 

    """
    Convert values in xarray DataArray from hPa to Pa

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """
    
    # do conversion
    da = da * 100 

    # update units attribute
    da["units"] = "Pa" 

    return da

def _convert_Pa_to_hPa(da: xr.DataArray): 

    """
    Convert values in xarray DataArray from Pa to hPa

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """
    
    # do conversion
    da = da / 100 

    # update units attribute
    da["units"] = "hPa" 

    return da


# better function name welcome! 
def _convert_mm_to_kg_m2_s1(da: xr.DataArray): 

    """
    Convert daily (!) values in xarray DataArray from mm to kg m^-2 s^-1

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray to convert

    Returns
    -------
    xr.DataArray
        The  converted xarray DataArray
    """
    
    # first, get timestep (frequency) by calculating the difference between the first consecutive time values in seconds
    timestep_nseconds = da.time.diff(dim='time').values[0]/  np.timedelta64(1, 's')

    # do conversion
    da = da / timestep_nseconds # mm to kg m^-2 s^-1  

    # update units attribute
    da["units"] = "kg m-2 s-1" 

    return da
