from valenspy._utilities import load_yml
from valenspy._unit_conversions import UNIT_CONVERSION_FUNCTIONS, EQUIVALENT_UNITS

import pytest

@pytest.mark.parametrize("lookup_name", ["EOBS_lookup", "CCLM_lookup", "ALARO-SFX_K_lookup","ERA5_lookup"])
def test_lookup_unit_conversion_coverage(lookup_name):
    """
    Test to see if there are any units in the lookup table which should be converted but are not defined in the UNIT_CONVERSION_FUNCTIONS or EQUIVALENT_UNITS.
    """
    CORDEX_VARIABLES = load_yml("CORDEX_variables")
    cf_units = set([var_attr["units"] for _, var_attr in CORDEX_VARIABLES.items()])

    equivalent_units = set([unit for unit in EQUIVALENT_UNITS])
    convertable_units = set([unit for unit in UNIT_CONVERSION_FUNCTIONS])
    all_convertable_units = equivalent_units.union(convertable_units)

    lookup_table = load_yml(lookup_name)
    to_convert_raw_units = set([var_attr["raw_units"] for  var, var_attr in lookup_table.items() if var_attr["raw_units"] != CORDEX_VARIABLES[var]["units"]])
    
    assert to_convert_raw_units.issubset(all_convertable_units), f"{to_convert_raw_units - all_convertable_units} in {lookup_name}.yml should be converted for some variable but are not."

def test_overlap_between_CORDEX_units_and_convertable_units():
    """
    Test whether the units defined to be converted using the UNIT_CONVERSION_FUNCTIONS (or the EQUIVALENT_UNITS) have no overlap with the units defined in the CORDEX variables.
    This avoids incorectly converting units which are already in the CORDEX variables.
    """
    equivalent_units = set([unit for unit in EQUIVALENT_UNITS])
    convertable_units = set([unit for unit in UNIT_CONVERSION_FUNCTIONS])
    all_convertable_units = equivalent_units.union(convertable_units)

    CORDEX_VARIABLES = load_yml("CORDEX_variables")
    cf_units = set([var_attr["units"] for _, var_attr in CORDEX_VARIABLES.items()])
    
    assert cf_units.isdisjoint(all_convertable_units), f"{cf_units.intersection(all_convertable_units)} are in the CORDEX_VARIABLE units and the equivalent_units or UNIT_CONVERSION_FUNCTIONS. This will lead to incorrect unit conversion."

