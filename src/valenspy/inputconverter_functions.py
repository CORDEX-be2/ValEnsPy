from pathlib import Path
import xarray as xr
from yaml import safe_load
from valenspy.cf_checks import is_cf_compliant, cf_status

# get path of source code (current path)
src_path = Path(__file__).resolve().parent

# open CORDEX variable lookup dictionary
with open(src_path / "ancilliary_data" / "CORDEX_variables.yml") as file:
    CORDEX_VARIABLES = safe_load(file)


def _non_convertor(paths,  metadata_info=None):
    """A dummy function that does not convert the input paths."""
    return xr.open_mfdataset(paths, combine="by_coords", chunks="auto")


def EOBS_to_CF(paths,  metadata_info=None) -> xr.Dataset:
    """
    Convert the EOBS netCDF paths to an xarray netCDF in CF convention

    Parameters
    ----------
    paths : Path or a list of Paths
        The path(s) to the netCDF file of specific variable to convert

    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """

    # open the observation dataset
    ds = xr.open_mfdataset(paths, combine="by_coords", chunks="auto")

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator?
    obsdata_name = "EOBS"

    with open(
        src_path / "ancilliary_data" / Path(obsdata_name + "_lookup.yml")
    ) as file:
        obs_LOOKUP = safe_load(file)

    # make EOBS CF compliant

    for var_obs in ds.data_vars:

        # Get the CORDEX variable in the observational dataset using the observational lookup table
        var = next(
            (k for k, v in obs_LOOKUP.items() if v.get("obs_name") == var_obs), None
        )

        if var:  # Dont processes variables that are not in the lookup table.

            # update variable name to CORDEX variable name
            ds = ds.rename_vars({obs_LOOKUP[var]["obs_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded EOBS units for units different to CORDEX
            if obs_LOOKUP[var]["obs_units"] == "Celcius":
                ds[var] = _convert_Celcius_to_Kelvin(ds[var])
            elif obs_LOOKUP[var]["obs_units"] == "hPa":
                ds[var] = _convert_hPa_to_Pa(ds[var])  # hPa to Pa
            elif (
                obs_LOOKUP[var]["obs_units"] == "mm"
            ):  # ! note observations remain daily time frequency
                ds[var] = _convert_mm_to_kg_m2s(ds[var])  # mm to kg m^-2 s^-1

            # add necessary metadata
            ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var][
                "standard_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["long_name"] = CORDEX_VARIABLES[var][
                "long_name"
            ]  # from the CORDEX look-up table

            ds[var].attrs["original_name"] = obs_LOOKUP[var]["obs_name"]
            ds[var].attrs["original_long_name"] = obs_LOOKUP[var]["obs_long_name"]

            # rename dimensions if not yet renamed
            if "lon" not in ds.coords: 
                ds = ds.rename({"longitude": "lon"})
            if "lat" not in ds.coords: 
                ds = ds.rename({"latitude": "lat"})

            # convert the time dimension to a pandas datetime index
            ds[var]["time"] = pd.to_datetime(ds[var].time)

            # additional attributes, on data-array level -- hard coded for EOBS
            ds[var].attrs["dataset"] = obsdata_name
            
            # if metadata_info is given, create global attributes
            if metadata_info: 
                for key, value in metadata_info.items(): 
                    ds[var].attrs[key] = value

            # if not, include hard-coded attributes (dataset dependent!)
            else: 
                ds[var].attrs["freq"] = "daily"
                ds[var].attrs["spatial_resolution"] = "0.1deg"
                ds[var].attrs["region"] = "europe"
    


    # set global attributes for whole dataset
    ds.attrs["dataset"] = obsdata_name

    # if metadata_info is given, create global attributes
    if metadata_info: 
        for key, value in metadata_info.items(): 
            ds.attrs[key] = value

    # if not, include hard-coded attributes (dataset dependent!)
    else: 
        ds.attrs["freq"] = "daily"  
        ds.attrs["spatial_resolution"] = "0.1deg"
        ds.attrs["region"] = "europe"


    # Soft check for CF compliance
    cf_status(ds)

    return ds


def ERA5_to_CF(file: Path, metadata_info=None) -> Path:
    """
    Convert the ERA5 netCDF file to a xarray Dataset in CF convention

    Parameters
    ----------
    file : Path
        The path to the netCDF file of specific variable to convert

    Returns
    -------
    Dataset
        The CF compliant ERA5 observations for the specified variable.
    """

    obsdata_name = "ERA5"


    # open the observation dataset
    ds = xr.open_mfdataset(file, combine="by_coords", chunks="auto")

    # open observational specific lookup dictionary
    with open(src_path / "ancilliary_data" / Path("ERA5_lookup.yml")) as lookup_file:
        obs_LOOKUP = safe_load(lookup_file)

    # make observation CF compliant
    for var_obs in ds.data_vars:

        # Get the CORDEX variable in the observational dataset using the observational lookup table
        var = next(
            (k for k, v in obs_LOOKUP.items() if v.get("obs_name") == var_obs), None
        )

        if var:  # Dont processes variables that are not in the lookup table.

            # update variable name to CORDEX variable name
            ds = ds.rename_vars({obs_LOOKUP[var]["obs_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded ERA5 units for CORDEX CORE, double check beyond.
            if (obs_LOOKUP[var]["obs_units"] == "Celcius") or (
                obs_LOOKUP[var]["obs_units"] == "degC"
            ):
                ds[var] = _convert_Celcius_to_Kelvin(ds[var])

            elif obs_LOOKUP[var]["obs_units"] == "hPa":
                ds[var] = _convert_hPa_to_Pa(ds[var])  # hPa to Pa

            elif (obs_LOOKUP[var]["obs_units"] == "mm") or (
                obs_LOOKUP[var]["obs_units"] == "mm/hr"
            ):
                ds[var] = _convert_mm_to_kg_m2s(
                    ds[var]
                )  # mm to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif (obs_LOOKUP[var]["obs_units"] == "m") or (
                obs_LOOKUP[var]["obs_units"] == "m/hr"
            ):
                ds[var] = _convert_m_to_kg_m2s(
                    ds[var]
                )  # m to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif obs_LOOKUP[var]["obs_units"] == "J/m^2":
                ds[var] = _convert_m_to_kg_m2s(
                    ds[var]
                )  # m to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion_convert_J_m2_to_W_m2

            # add necessary metadata
            ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var][
                "standard_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["long_name"] = CORDEX_VARIABLES[var][
                "long_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["original_name"] = obs_LOOKUP[var]["obs_name"]
            ds[var].attrs["original_long_name"] = obs_LOOKUP[var]["obs_long_name"]

            # rename dimensions if not yet renamed
            if "lon" not in ds[var].coords: 
                ds = ds.rename({"longitude": "lon"})
            if "lat" not in ds[var].coords: 
                ds = ds.rename({"latitude": "lat"})

            # convert the time dimension to a pandas datetime index --  do we want this to happen within the convertor? Or do we leave it up to the user?
            ds[var]["time"] = pd.to_datetime(ds[var].time)

            # additional attributes at data array level.
            ds[var].attrs["dataset"] = obsdata_name

            if metadata_info: 
                for key, value in metadata_info.items(): 
                    ds[var].attrs[key] = value

            # if not, include hard-coded attributes (dataset dependent!)
            else: 
                ds[var].attrs["freq"] = _determine_time_interval(ds[var])
    


    # set attributes in whole dataset 
    ds.attrs["dataset"] = obsdata_name

    # if metadata_info is given, create global attributes
    if metadata_info: 
        for key, value in metadata_info.items(): 
            ds.attrs[key] = value

    # if not, include hard-coded attributes (dataset dependent!)
    else: 
        ds.attrs["freq"] = _determine_time_interval(ds) # automatically check on time interval based on time axis. 

    # Soft check for CF compliance
    cf_status(ds)

    return ds



def ERA5Land_to_CF(file: Path,  metadata_info=None) -> Path:
    """
    Convert the ERA5-Land netCDF file to a xarray Dataset in CF convention

    Parameters
    ----------
    file : Path
        The path to the netCDF file of specific variable to convert

    Returns
    -------
    Dataset
        The CF compliant ERA5-Land observations for the specified variable.
    """

    obsdata_name = "ERA5-Land"

    # open the observation dataset
    ds = xr.open_mfdataset(file, combine="by_coords", chunks="auto")

    # open observational specific lookup dictionary
    with open(src_path / "ancilliary_data" / Path("ERA5_lookup.yml")) as lookup_file:
        obs_LOOKUP = safe_load(lookup_file)

    # make observation CF compliant
    for var_obs in ds.data_vars:

        # Get the CORDEX variable in the observational dataset using the observational lookup table
        var = next(
            (k for k, v in obs_LOOKUP.items() if v.get("obs_name") == var_obs), None
        )

        if var:  # Dont processes variables that are not in the lookup table.

            # update variable name to CORDEX variable name
            ds = ds.rename_vars({obs_LOOKUP[var]["obs_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded ERA5 units for CORDEX CORE, double check beyond.
            if (obs_LOOKUP[var]["obs_units"] == "Celcius") or (
                obs_LOOKUP[var]["obs_units"] == "degC"
            ):
                ds[var] = _convert_Celcius_to_Kelvin(ds[var])

            elif obs_LOOKUP[var]["obs_units"] == "hPa":
                ds[var] = _convert_hPa_to_Pa(ds[var])  # hPa to Pa

            elif (obs_LOOKUP[var]["obs_units"] == "mm") or (
                obs_LOOKUP[var]["obs_units"] == "mm/hr"
            ):
                ds[var] = _convert_mm_to_kg_m2s(
                    ds[var]
                )  # mm to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif (obs_LOOKUP[var]["obs_units"] == "m") or (
                obs_LOOKUP[var]["obs_units"] == "m/hr"
            ):
                ds[var] = _convert_m_to_kg_m2s(
                    ds[var]
                )  # m to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif obs_LOOKUP[var]["obs_units"] == "J/m^2":
                ds[var] = _convert_J_m2_to_W_m2(
                    ds[var]
                )  # J to W m2 conversion function reads time frequency (nseconds) of input ds to do conversion_convert_J_m2_to_W_m2

            # add necessary metadata
            ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var][
                "standard_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["long_name"] = CORDEX_VARIABLES[var][
                "long_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["original_name"] = obs_LOOKUP[var]["obs_name"]
            ds[var].attrs["original_long_name"] = obs_LOOKUP[var]["obs_long_name"]

            # rename dimensions if not yet renamed
            if "lon" not in ds.coords: 
                ds = ds.rename({"longitude": "lon"})
            if "lat" not in ds.coords: 
                ds = ds.rename({"latitude": "lat"})

            # convert the time dimension to a pandas datetime index --  do we want this to happen within the convertor? Or do we leave it up to the user?
            ds[var]["time"] = pd.to_datetime(ds[var].time)

            # additional attributes -- set both globally at dataset level as at data array level
            ds[var].attrs["dataset"] = obsdata_name

            if metadata_info: 
                for key, value in metadata_info.items(): 
                    ds[var].attrs[key] = value

            # if not, include hard-coded attributes (dataset dependent!)
            else: 
                ds[var].attrs["freq"] = _determine_time_interval(ds[var])

    # set attributes in whole dataset 
    ds.attrs["dataset"] = obsdata_name

    # if metadata_info is given, create global attributes
    if metadata_info: 
        for key, value in metadata_info.items(): 
            ds.attrs[key] = value

    # if not, include hard-coded attributes (dataset dependent!)
    else: 
        ds.attrs["freq"] = _determine_time_interval(ds) # automatically check on time interval based on time axis. 

    # Soft check for CF compliance
    cf_status(ds)

    return ds


def CLIMATE_GRID_to_CF(file: Path,  metadata_info=None) -> xr.Dataset:
    """
    Convert the CLIMATE_GRID netCDF paths to an xarray netCDF in CF convention

    Parameters
    ----------
    paths : Path or a list of Paths
        The path(s) to the netCDF file of specific variable to convert

    Returns
    -------
    Dataset
        The CF compliant CLIMATE_GRID observations for the specified variable.
    """

    obsdata_name = "CLIMATE_GRID"       
    
    # open the observation dataset
    ds = xr.open_mfdataset(file, combine="by_coords", chunks="auto")

    # open observational specific lookup dictionary
    with open(src_path / "ancilliary_data" / Path(obsdata_name+"_lookup.yml")) as lookup_file:
        obs_LOOKUP = safe_load(lookup_file)

    # make observation CF compliant
    for var_obs in ds.data_vars:

        # Get the CORDEX variable in the observational dataset using the observational lookup table
        var = next(
            (k for k, v in obs_LOOKUP.items() if v.get("obs_name") == var_obs), None
        )

        if var:  # Dont processes variables that are not in the lookup table.
            # update variable name to CORDEX variable name
            ds = ds.rename_vars({obs_LOOKUP[var]["obs_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded ERA5 units for CORDEX CORE, double check beyond.
            if (obs_LOOKUP[var]["obs_units"] == "Celcius") or (
                obs_LOOKUP[var]["obs_units"] == "degC"
            ):
                ds[var] = _convert_Celcius_to_Kelvin(ds[var])

            elif obs_LOOKUP[var]["obs_units"] == "hPa":
                ds[var] = _convert_hPa_to_Pa(ds[var])  # hPa to Pa

            elif (obs_LOOKUP[var]["obs_units"] == "mm") or (
                obs_LOOKUP[var]["obs_units"] == "mm/hr"
            ):
                ds[var] = _convert_mm_to_kg_m2s(
                    ds[var]
                )  # mm to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif (obs_LOOKUP[var]["obs_units"] == "m") or (
                obs_LOOKUP[var]["obs_units"] == "m/hr"
            ):
                ds[var] = _convert_m_to_kg_m2s(
                    ds[var]
                )  # m to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif obs_LOOKUP[var]["obs_units"] == "kWh/m2/day":
                ds[var] = _convert_kWh_m2_day_to_W_m2(
                    ds[var]
                ) 

            # add necessary metadata
            ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var][
                "standard_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["long_name"] = CORDEX_VARIABLES[var][
                "long_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["original_name"] = obs_LOOKUP[var]["obs_name"]
            ds[var].attrs["original_long_name"] = obs_LOOKUP[var]["obs_long_name"]

            # regrid to lat lon 

            # convert the time dimension to a pandas datetime index --  do we want this to happen within the convertor? Or do we leave it up to the user?
            ds[var]["time"] = pd.to_datetime(ds[var].time)    
            # additional attributes -- set both globally at dataset level as at data array level
            ds[var].attrs["dataset"] = obsdata_name

            if metadata_info: 
                for key, value in metadata_info.items(): 
                    ds[var].attrs[key] = value

            # if not, include hard-coded attributes (dataset dependent!)
            else: 
                ds[var].attrs["freq"] = "daily"
                ds[var].attrs["spatial_resolution"] = "~5km, regridded to regular lat, lon grid"
                ds[var].attrs["region"] = "belgium" # leave empty per default. 
    


    # set attributes in whole dataset 
    ds.attrs["dataset"] = obsdata_name

    # if metadata_info is given, create global attributes
    if metadata_info: 
        for key, value in metadata_info.items(): 
            ds[var].attrs[key] = value

    # if not, include hard-coded attributes (dataset dependent!)
    else: 
        ds.attrs["freq"] = "daily"
        ds.attrs["spatial_resolution"] = "5 km"
        ds.attrs["region"] = "belgium" # leave empty per default. 

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
    da = da + 273.15  # Celcius to Kelvin

    # update units attribute --  naming of units defined in ./src/valenspy/ancilliary_data/CORDEX_variables.yml
    da.attrs["units"] = "K"

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
    da = da - 273.15  # Kelvin to Celcius

    # update units attribute
    da.attrs["units"] = "°C"

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
    da.attrs["units"] = "Pa"

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
    da.attrs["units"] = "hPa"

    return da


def _convert_mm_to_kg_m2s(da: xr.DataArray):
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
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    # do conversion
    da = da / timestep_nseconds  # mm to kg m^-2 s^-1

    # update units attribute
    da.attrs["units"] = "kg m-2 s-1"

    return da


def _convert_m_to_kg_m2s(da: xr.DataArray):
    """
    Convert values in xarray DataArray from mm hr^-1 to kg m^-2 s^-1

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
    da = da * 1000 / 3600  # mm hr^-1 to kg m^-2 s^-1

    # update units attribute
    da.attrs["units"] = "kg m-2 s-1"

    return da

    # Convert J/m²/hr to W/m²
    w_m2 = j_m2_hr / seconds_per_hour


def _convert_J_m2_to_W_m2(da: xr.DataArray):
    """
    Convert values in xarray DataArray from J m^2 to W m^2

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
    timestep_nseconds = da.time.diff(dim="time").values[0] / np.timedelta64(1, "s")

    # do conversion
    da = da / timestep_nseconds  # J m^2 to W m^2

    # update units attribute
    da.attrs["units"] = "W m-2"

    return da

def _convert_kWh_m2_day_to_W_m2(da: xr.DataArray):
    """
    Convert values in xarray DataArray from kWh/m2/day to W m^2

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
    da = da *1000 / 86400   # kWh/m2/day to W m^2

    # update units attribute
    da.attrs["units"] = "W m-2"

    return da

def _determine_time_interval(da: xr.DataArray): 

    """
    Find the time interval (freq) of the input data array based on it's time axis, by calculating the difference between the first two time instances. 

    Parameters
    ----------
    da : xr.DataArray
        The xarray DataArray with time axis to check the time interval

    Returns
    -------
    freq : string
        The frequency string containing "hourly, daily, monthly or yearly"
    """

    diff = da.time.diff(dim="time").values[0]

    # Check for exact differences
    if diff == np.timedelta64(1, 'h'):
        freq="hourly"
    elif diff == np.timedelta64(1, 'D'):
        freq="daily"
    elif diff == np.timedelta64(1, 'M'):
        freq="monthly"
    elif diff == np.timedelta64(1, 'Y'):
        freq="yearly"
    else:
        return "Difference does not match exact hourly, daily, monthly, or yearly units."

    return freq
