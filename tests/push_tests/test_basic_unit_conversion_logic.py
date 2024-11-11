from valenspy._utilities import load_yml
from valenspy.input.unit_converter import UNIT_CONVERSION_FUNCTIONS, EQUIVALENT_UNITS

import pytest

lookup_tables = ["EOBS_lookup", "CCLM_lookup", "ALARO-SFX_K_lookup","ERA5_lookup","CLIMATE_GRID_lookup", "RADCLIM_lookup", "GPM_IMERG_lookup"]

@pytest.mark.parametrize("lookup_name", lookup_tables)
def test_lookup_unit_conversion_coverage(lookup_name):
    """
    Test to see if there are any units in the lookup table which should be converted but are not defined in the UNIT_CONVERSION_FUNCTIONS or EQUIVALENT_UNITS.
    The UNIT_CONVERSION_FUNCTIONS contain the mapping of units to the conversion functions - indicating that these units can be converted.
    The EQUIVALENT_UNITS contain the mapping of units that are equivalent to each other - indicating that these units can be converted to equivalent units which are defined in the UNIT_CONVERSION_FUNCTIONS.
    """
    CORDEX_VARIABLES = load_yml("CORDEX_variables")
    cf_units = set([var_attr["units"] for _, var_attr in CORDEX_VARIABLES.items()])

    equivalent_units = set([unit for unit in EQUIVALENT_UNITS])
    convertable_units = set([unit for unit in UNIT_CONVERSION_FUNCTIONS])
    all_convertable_units = equivalent_units.union(convertable_units)

    lookup_table = load_yml(lookup_name)
    to_convert_raw_units = set([var_attr["raw_units"] for  var, var_attr in lookup_table.items() if var_attr["raw_units"] != CORDEX_VARIABLES[var]["units"]])
    
    assert to_convert_raw_units.issubset(all_convertable_units), f"{to_convert_raw_units - all_convertable_units} in {lookup_name}.yml should be converted for some variable but are not."

@pytest.mark.parametrize("lookup_name", lookup_tables)
def test_lookup_raw_name_uniqueness(lookup_name):
    """
    Test whether the raw names in each lookup tables are unique, either as a single string or a unique combination of strings and aggregation functions.
    """
    lookup_table = load_yml(lookup_name)

    raw_names = [var_attr["raw_name"] for var_attr in lookup_table.values()]
    unique_raw_names = set(raw_names)
    non_unique_raw_names = set([raw_name for raw_name in raw_names if raw_names.count(raw_name) > 1])

    assert len(non_unique_raw_names) == 0, f"Raw names {non_unique_raw_names} in {lookup_name}.yml are not unique."

def test_equivalent_units_definition():
    """
    Test whether the equivalent units are defined correctly in the EQUIVALENT_UNITS dictionary.
    The keys should be units that are not defined in the UNIT_CONVERSION_FUNCTIONS or the CORDEX_VARIABLES_units.
    The values should be the units that are defined in the CORDEX_VARIABLES or the UNIT_CONVERSION_FUNCTIONS.
    """
    equivalent_units = set([unit for unit in EQUIVALENT_UNITS])
    convertable_units = set([unit for unit in UNIT_CONVERSION_FUNCTIONS])
    CORDEX_VARIABLES_units = set(
        [var_attr["units"] for _, var_attr in load_yml("CORDEX_variables").items()]
    )

    convertable_or_CORDEX_units = convertable_units.union(CORDEX_VARIABLES_units)

    for key, value in EQUIVALENT_UNITS.items():
        assert key not in convertable_or_CORDEX_units, f"Key: {key} in EQUIVALENT_UNITS should not be defined in UNIT_CONVERSION_FUNCTIONS or CORDEX_VARIABLES_units."
        assert value in convertable_or_CORDEX_units, f"Value: {value} in EQUIVALENT_UNITS should be defined in UNIT_CONVERSION_FUNCTIONS or CORDEX_VARIABLES_units."