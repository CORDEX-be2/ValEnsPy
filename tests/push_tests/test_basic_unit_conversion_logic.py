from valenspy._utilities import load_yml

import xclim
import pytest

lookup_tables = ["EOBS_lookup", "CCLM_lookup", "ALARO-SFX_K_lookup","ERA5_lookup","CLIMATE_GRID_lookup", "RADCLIM_lookup"]

@pytest.mark.parametrize("lookup_name", lookup_tables)
def test_lookup_unit_conversion(lookup_name):
    """
    Test whether the units in the lookup tables can be converted using xclim's units module.
    """
    CORDEX_VARIABLES = load_yml("CORDEX_variables")
    lookup_table = load_yml(lookup_name)
    
    for var, var_attr in lookup_table.items():
        raw_units = var_attr["raw_units"]
        units = CORDEX_VARIABLES[var].get("units")
        try:
            xclim.core.units.convert_units_to(raw_units, units, context="infer")
        except Exception as e:
            pytest.fail(f"For {var}: Units {raw_units} in {lookup_name}.yml for variable {var} cannot be converted to {units}. Error: {e}")