from pathlib import Path
from valenspy.cf_checks import is_cf_compliant, cf_status
from valenspy._utilities import load_yml
from valenspy._unit_conversions import convert_all_units_to_CF, _convert_Celcius_to_Kelvin, _convert_hPa_to_Pa, _convert_mm_to_kg_m2s, _convert_m_to_kg_m2s, _convert_J_m2_to_W_m2
import xarray as xr
import pandas as pd
import numpy as np

CORDEX_VARIABLES = load_yml("CORDEX_variables")

def _set_global_attributes(ds: xr.Dataset, metadata_info):
    for key, value in metadata_info.items():
        ds.attrs[key] = value
    
    return ds

def EOBS_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the EOBS xarray netCDF to a CF compliant xarray Dataset

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of EOBS observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """
    #Unique information regarding the EOBS dataset
    obsdata_name = "EOBS"
    raw_LOOKUP = load_yml(f"{obsdata_name}_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq":"daily", "spatial_resolution":"0.1deg", "region":"europe"}

    #Convert all units to CF, add metadata and set global attributes
    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)

    # Soft check for CF compliance
    cf_status(ds)

    return ds

def ERA5_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the ERA5 xarray dataset to a xarray Dataset in CF convention

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ERA5 observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant ERA5 observations for the specified variable.
    """
    #Unique information regarding the ERA5 dataset
    obsdata_name = "ERA5"
    raw_LOOKUP = load_yml(f"{obsdata_name}_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq": _determine_time_interval(ds)}

    #Convert all units to CF, add metadata and set global attributes
    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)
    
    cf_status(ds)

    return ds

def ERA5Land_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the ERA5-Land xarray dataset to a xarray Dataset in CF convention.
    Uses the same lookup table as ERA5.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ERA5-Land observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant ERA5-Land observations for the specified variable.
    """
    #Unique information regarding the ERA5 dataset
    obsdata_name = "ERA5-Land"
    raw_LOOKUP = load_yml(f"ERA5_lookup") #Note that the ERA5 lookup is used for both ERA5 and ERA5-Land

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq": _determine_time_interval(ds)}

    #Convert all units to CF, add metadata and set global attributes
    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)
    
    cf_status(ds)

    return ds

def CLIMATE_GRID_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the CLIMATE_GRID xarray dataset to a xarray Dataset in CF convention

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CLIMATE_GRID observations to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant CLIMATE_GRID observations for the specified variable.
    """
    obsdata_name = "CLIMATE_GRID"
    raw_LOOKUP = load_yml(f"{obsdata_name}_lookup")

    if metadata_info is None: #Set standard metadata if not provided
        metadata_info = {"freq":"daily", "spatial_resolution":"0.07° x 0.045° (~5km)", "region":"belgium"}

    metadata_info["dataset"] = obsdata_name

    ds = convert_all_units_to_CF(ds, raw_LOOKUP, metadata_info)
    ds = _set_global_attributes(ds, metadata_info)
    
    cf_status(ds)

    return ds

def CCLM_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the CCLM xarray netCDF to a CF compliant xarray Dataset

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CCLM simulation to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator?
    model_name = "CCLM"
    src_path = Path("../src/valenspy")

    # open observational specific lookup dictionary
    raw_LOOKUP = load_yml(model_name + "_lookup")

    # make model dataset CF compliant
    for var_mod in ds.data_vars:

        # Get the CORDEX variable in the model dataset using the model-specific lookup table
        var = next(
            (k for k, v in raw_LOOKUP.items() if v is not None and v.get("raw_name") == var_mod),
            None
        )
        if var:  # Dont processes variables that are not in the lookup table.

            # update variable name to CORDEX variable name
            ds = ds.rename_vars({raw_LOOKUP[var]["raw_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded ERA5 units for CORDEX CORE, double check beyond.
            if (raw_LOOKUP[var]["raw_units"] == "Celcius") or (
                raw_LOOKUP[var]["raw_units"] == "degC"
            ):
                ds[var] = _convert_Celcius_to_Kelvin(ds[var])

            elif raw_LOOKUP[var]["raw_units"] == "hPa":
                ds[var] = _convert_hPa_to_Pa(ds[var])  # hPa to Pa

            elif (raw_LOOKUP[var]["raw_units"] == "mm") or (
                raw_LOOKUP[var]["raw_units"] == "mm/hr"
            ):
                ds[var] = _convert_mm_to_kg_m2s(
                    ds[var]
                )  # mm to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif raw_LOOKUP[var]["raw_units"] == "kWh/m2/day":
                ds[var] = _convert_kWh_m2_day_to_W_m2(
                    ds[var]
                )  # kWh/m2/day to W m^-2 conversion function reads time frequency (nseconds) of input ds to do conversion_convert_J_m2_to_W_m2

            elif raw_LOOKUP[var]["raw_units"] == "m/s":
                ds[var].attrs["units"] = CORDEX_VARIABLES[var][
                    "units"
                ]  # put units as m s-1

            # add necessary metadata
            ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var][
                "standard_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["long_name"] = CORDEX_VARIABLES[var][
                "long_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["original_name"] = raw_LOOKUP[var]["raw_name"]
            ds[var].attrs["original_long_name"] = raw_LOOKUP[var]["raw_long_name"]

            # convert the time dimension to a pandas datetime index
            ds[var]["time"] = pd.to_datetime(ds[var].time)

            # additional attributes -- set both globally at dataset level as at data array level
            ds[var].attrs["dataset"] = model_name

            if metadata_info:
                for key, value in metadata_info.items():
                    ds[var].attrs[key] = value

 

    # set attributes in whole dataset
    ds.attrs["dataset"] = model_name

    # if metadata_info is given, create global attributes
    if metadata_info:
        for key, value in metadata_info.items():
            ds[var].attrs[key] = value

    # Soft check for CF compliance
    cf_status(ds)

    return ds

def ALARO_K_to_CF(ds: xr.Dataset, metadata_info=None) -> xr.Dataset:
    """
    Convert the ALARO_K (converted with Kwintens R scripts) from xarray netCDF to a CF compliant xarray Dataset

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ALARO_K simulation to convert
    metadata_info : dict, optional
        A dictionary containing additional metadata information to add to the dataset

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator?
    model_name = "ALARO-SFX_K"

    # open observational specific lookup dictionary
    raw_LOOKUP = load_yml(model_name + "_lookup")

    # make model dataset CF compliant
    for var_mod in ds.data_vars:

        # Get the CORDEX variable in the model dataset using the model-specific lookup table
        var = next(
            (k for k, v in raw_LOOKUP.items() if v is not None and v.get("raw_name") == var_mod),
            None
        )
        if var:  # Dont processes variables that are not in the lookup table.

            # update variable name to CORDEX variable name
            ds = ds.rename_vars({raw_LOOKUP[var]["raw_name"]: var})

            # from here on, use CORDEX variable name to access data array and do rest of conversion

            # Unit conversion - hard coded ERA5 units for CORDEX CORE, double check beyond.
            if (raw_LOOKUP[var]["raw_units"] == "Celcius") or (
                raw_LOOKUP[var]["raw_units"] == "degC"
            ):
                ds[var] = _convert_Celcius_to_Kelvin(ds[var])

            elif raw_LOOKUP[var]["raw_units"] == "hPa":
                ds[var] = _convert_hPa_to_Pa(ds[var])  # hPa to Pa

            elif (raw_LOOKUP[var]["raw_units"] == "mm") or (
                raw_LOOKUP[var]["raw_units"] == "mm/hr"
            ):
                ds[var] = _convert_mm_to_kg_m2s(
                    ds[var]
                )  # mm to kg m^-2 s^-1 conversion function reads time frequency (nseconds) of input ds to do conversion

            elif raw_LOOKUP[var]["raw_units"] == "kWh/m2/day":
                ds[var] = _convert_kWh_m2_day_to_W_m2(
                    ds[var]
                )  # kWh/m2/day to W m^-2 conversion function reads time frequency (nseconds) of input ds to do conversion_convert_J_m2_to_W_m2

            elif raw_LOOKUP[var]["raw_units"] == "m/s":
                ds[var].attrs["units"] = CORDEX_VARIABLES[var][
                    "units"
                ]  # put units as m s-1

            # add necessary metadata
            ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var][
                "standard_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["long_name"] = CORDEX_VARIABLES[var][
                "long_name"
            ]  # from the CORDEX look-up table
            ds[var].attrs["original_name"] = raw_LOOKUP[var]["raw_name"]
            ds[var].attrs["original_long_name"] = raw_LOOKUP[var]["raw_long_name"]

            # convert the time dimension to a pandas datetime index
            ds[var]["time"] = pd.to_datetime(ds[var].time)

            # additional attributes -- set both globally at dataset level as at data array level
            ds[var].attrs["dataset"] = model_name

            if metadata_info:
                for key, value in metadata_info.items():
                    ds[var].attrs[key] = value

    if "rain_convective" in ds.data_vars and "rain_stratiform" in ds.data_vars:
        ds["pr"] = _convert_mm_to_kg_m2s(ds["rain_convective"] + ds["rain_stratiform"])
        ds["pr"].attrs["standard_name"] = "precipitation_flux"
        ds["pr"].attrs["long_name"] = "Precipitation"
        ds["pr"].attrs["dataset"] = model_name
        ds["pr"].attrs["original_name"] = "rain_convective + rain_stratiform"
        if metadata_info:
                for key, value in metadata_info.items():
                    ds[var].attrs[key] = value

    # set attributes in whole dataset
    ds.attrs["dataset"] = model_name

    # if metadata_info is given, create global attributes
    if metadata_info:
        for key, value in metadata_info.items():
            ds[var].attrs[key] = value

    # Soft check for CF compliance
    cf_status(ds)

    return ds

# helper functions - can be moved to more appropriate place

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
    if diff == np.timedelta64(1, "h"):
        freq = "hourly"
    elif diff == np.timedelta64(1, "D"):
        freq = "daily"
    elif diff == np.timedelta64(1, "M"):
        freq = "monthly"
    elif diff == np.timedelta64(1, "Y"):
        freq = "yearly"
    else:
        return (
            "Difference does not match exact hourly, daily, monthly, or yearly units."
        )

    return freq