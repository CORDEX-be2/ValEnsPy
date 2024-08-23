# helper functions for unit conversion - can be moved to more appropriate place

## imports for helpen functions
import xarray as xr
import pandas as pd
import warnings
from valenspy._utilities import load_yml
from valenspy.unit_conversion_functions import (
    _convert_Celcius_to_Kelvin,
    _convert_hPa_to_Pa,
    _convert_mm_to_kg_m2s,
    _convert_m_to_kg_m2s,
    _convert_J_m2_to_W_m2,
    _convert_kWh_m2_day_to_W_m2,
    _convert_fraction_to_percent,
    _determine_time_interval,
)

CORDEX_VARIABLES = load_yml("CORDEX_variables")

UNIT_CONVERSION_FUNCTIONS = {
    "Celcius": _convert_Celcius_to_Kelvin,
    "hPa": _convert_hPa_to_Pa,
    "mm": _convert_mm_to_kg_m2s,
    "mm/hr": _convert_m_to_kg_m2s,
    "m": _convert_m_to_kg_m2s,
    "m/hr": _convert_m_to_kg_m2s,
    "J/m^2": _convert_J_m2_to_W_m2,
    "kWh/m2/day": _convert_kWh_m2_day_to_W_m2,
    "1": _convert_fraction_to_percent
}

EQUIVALENT_UNITS = {
    "degC": "Celcius", 
    "m/s": "m s-1", 
    "(0 - 1)": "1"
    }

def convert_all_units_to_CF(ds: xr.Dataset, raw_LOOKUP, metadata_info: dict):
    """Convert all units for all variables in the dataset to the correct units by applying the correct conversion function.

    The raw_units attribute is used to determine the current units and therefore the conversion function to apply.
    The target units are defined in the CORDEX_variables.yml file. In addition it:
    - Handels variable metadata attributes specific to the var and conversion
    - Adds dataset metadata attributes to the variable attributes
    - Converts the time dimension to a pandas datetime object


    Parameters
    ----------
    ds : xr.Dataset
        The xarray dataset to convert
    raw_LOOKUP : dict
        The lookup table for the variables which corresponds the CORDEX variables to the variables in the dataset and their units
    metadata_info : dict
        The metadata information which should be added to the each variable in the dataset

    Returns
    -------
    xr.Dataset
        The converted xarray dataset
    """

    # Key: The unit of the raw data
    # Value: The unit of the CORDEX equivalent unit or the unit that is used to identify the conversion function


    for raw_var in ds.data_vars:
        var = next(
            (k for k, v in raw_LOOKUP.items() if v.get("raw_name") == raw_var), None
        )

        if var:  # Dont processes variables that are not in the lookup table.

            # rename variable to CORDEX variable name
            ds = ds.rename_vars({raw_var: var})

            # convert units - based on the raw units
            raw_units = raw_LOOKUP[var]["raw_units"]
            cordex_var_units = CORDEX_VARIABLES[var]["units"]

            # If the raw_units are in the equivalent_units, use the replacement unit
            if raw_units in EQUIVALENT_UNITS:
                raw_units = EQUIVALENT_UNITS[raw_units]

            if raw_units in UNIT_CONVERSION_FUNCTIONS:

                ds[var] = UNIT_CONVERSION_FUNCTIONS[raw_units](
                    ds[var]
                )  # Do the conversion
            elif not raw_units == cordex_var_units:
                # Throw a warning that the raw_unit in the lookup table is not implemented
                warnings.warn(
                    f"Unit conversion for {raw_units} to {cordex_var_units} is not implemented for variable {var}."
                )

            if var in ds:  # If the renaming occured add the metadata attributes
                # Metadata attributes
                ds[var].attrs["standard_name"] = CORDEX_VARIABLES[var]["standard_name"]
                ds[var].attrs["long_name"] = CORDEX_VARIABLES[var]["long_name"]
                ds[var].attrs["units"] = CORDEX_VARIABLES[var]["units"]
                ds[var].attrs["original_name"] = raw_LOOKUP[var]["raw_name"]
                ds[var].attrs["original_long_name"] = raw_LOOKUP[var]["raw_long_name"]
                ds[var].attrs["original_units"] = raw_LOOKUP[var]["raw_units"]

                ds[var]["time"] = pd.to_datetime(ds[var].time)

                if metadata_info:
                    for key, value in metadata_info.items():
                        ds[var].attrs[key] = value

                # If freq is not in the metadata_info, we can try to infer it from the time dimension
                if "freq" not in ds[var].attrs:
                    ds[var].attrs["freq"] = _determine_time_interval(ds[var])

    return ds
