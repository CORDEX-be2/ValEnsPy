from valenspy._utilities import load_yml
from valenspy._unit_conversions import unit_conversion_functions, EQUIVALENT_UNITS

import pytest

@pytest.mark.parametrize("lookup_name", ["EOBS_lookup", "CCLM_lookup", "ALARO-SFX_K_lookup","ERA5_lookup"])
def test_lookup_unit_coverage(lookup_name):
    """
    Testing whether the units defined in a lookup table are either:
    1. Already in one of the CF defined conversion units
    2. Will be converted to one of the CF defined conversion units as it is defined in
        a. The unit_conversion_functions - a dictionary linking raw_units to the conversion function to the desired units
        b. The EQUIVALENT_UNITS - a dictionary listing all equivalent units which are linked to either a unit which can be converted to a CF unit or a CF unit itself
    """
    CORDEX_VARIABLES = load_yml("CORDEX_variables")

    cf_units = set([var_attr["units"] for _, var_attr in CORDEX_VARIABLES.items()])
    equivalent_units = set([unit for unit in EQUIVALENT_UNITS])
    convertable_units = set([unit for unit in unit_conversion_functions])
    covered_units = cf_units.union(equivalent_units).union(convertable_units)

    lookup_table = load_yml(lookup_name)
    raw_units = set([var_attr["raw_units"] for  _, var_attr in lookup_table.items()])
    assert raw_units.issubset(covered_units), f"The following units in {lookup_name}.yml are not covered: {raw_units - covered_units}"

def test_overlap_between_CORDEX_units_and_convertable_units():
    """
    Test whether the units defined to be converted using the unit_conversion_functions (or the EQUIVALENT_UNITS) have no overlap with the units defined in the CORDEX variables.
    This avoids incorectly converting units which are already in the CORDEX variables.
    """
    equivalent_units = set([unit for unit in EQUIVALENT_UNITS])
    convertable_units = set([unit for unit in unit_conversion_functions])
    all_convertable_units = equivalent_units.union(convertable_units)

    CORDEX_VARIABLES = load_yml("CORDEX_variables")
    cf_units = set([var_attr["units"] for _, var_attr in CORDEX_VARIABLES.items()])
    
    assert cf_units.isdisjoint(all_convertable_units), f"The following units are both in the CORDEX variables and the equivalent_units or unit_conversion_functions: {cf_units.intersection(all_convertable_units)}"