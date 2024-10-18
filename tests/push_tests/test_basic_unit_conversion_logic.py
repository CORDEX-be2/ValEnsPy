from valenspy._utilities import load_yml
from valenspy._unit_conversions import UNIT_CONVERSION_FUNCTIONS, EQUIVALENT_UNITS

import pytest


@pytest.mark.parametrize(
    "lookup_name",
    [
        "EOBS_lookup",
        "CCLM_lookup",
        "ALARO-SFX_K_lookup",
        "ERA5_lookup",
        "CLIMATE_GRID_lookup",
    ],
)
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
    to_convert_raw_units = set(
        [
            var_attr["raw_units"]
            for var, var_attr in lookup_table.items()
            if var_attr["raw_units"] != CORDEX_VARIABLES[var]["units"]
        ]
    )

    assert to_convert_raw_units.issubset(
        all_convertable_units
    ), f"{to_convert_raw_units - all_convertable_units} in {lookup_name}.yml should be converted for some variable but are not."


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
        assert (
            key not in convertable_or_CORDEX_units
        ), f"Key: {key} in EQUIVALENT_UNITS should not be defined in UNIT_CONVERSION_FUNCTIONS or CORDEX_VARIABLES_units."
        assert (
            value in convertable_or_CORDEX_units
        ), f"Value: {value} in EQUIVALENT_UNITS should be defined in UNIT_CONVERSION_FUNCTIONS or CORDEX_VARIABLES_units."
