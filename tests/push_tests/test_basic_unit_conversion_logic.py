from valenspy._utilities import load_yml

import xclim
import xarray as xr
import pandas as pd
import numpy as np
import pytest

lookup_tables = ["EOBS_lookup", "CCLM_lookup", "ALARO-SFX_K_lookup","ERA5_lookup","CLIMATE_GRID_lookup", "RADCLIM_lookup"]

@pytest.mark.parametrize("lookup_name", lookup_tables)
def test_lookup_unit_conversion(lookup_name):
    """
    Test whether the units in the lookup tables can be converted using xclim's units module.

    Test if conversion from raw_units (in the lookup table) to the units in CORDEX_variables.yml is possible.
    The raw_standard_name is used to infer the context for the conversions. 
    In particular see, https://xclim.readthedocs.io/en/stable/notebooks/units.html#Smart-conversions:-Precipitation 
    """
    CORDEX_VARIABLES = load_yml("CORDEX_variables")
    lookup_table = load_yml(lookup_name)
    
    for var, var_attr in lookup_table.items():
        raw_units = var_attr["raw_units"]
        units = CORDEX_VARIABLES[var].get("units")

        #Define the context for the conversion
        ctxs  = []
        ctxs.append(xclim.core.units.infer_context(CORDEX_VARIABLES[var].get("standard_name")))
        if var_attr.get("raw_standard_name"):
            ctxs.append(xclim.core.units.infer_context(var_attr.get("raw_standard_name")))

        context = "hydro" if "hydro" in ctxs else None
        
        if context == "hydro":
            # Smart conversions
            # https://xclim.readthedocs.io/en/stable/notebooks/units.html#Smart-conversions:-Precipitation 
            raw_standard_name = var_attr.get("raw_standard_name", None)
            da = xr.DataArray(
                np.random.rand(4, 1, 1), 
                dims=["time", "lat", "lon"], 
                coords={
                    "time": pd.date_range("2000-01-01", periods=4),
                    "lat": [0],
                    "lon": [0]
                },
                attrs={"units": raw_units, "standard_name": raw_standard_name}
            )
            try:
                xclim.core.units.convert_units_to(da, units)
            except Exception as e:
                pytest.fail(f"For {var}: Units {raw_units}, standard_name {raw_standard_name} in {lookup_name}.yml for variable {var} cannot be converted to {units}. Error: {e}")
        else:
            try:
                xclim.core.units.convert_units_to(raw_units, units)
            except Exception as e:
                pytest.fail(f"{var} - Units {raw_units} in {lookup_name}.yml cannot be converted to {units}. Error: {e}")