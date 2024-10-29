import valenspy as vp
from valenspy.cf_checks import is_cf_compliant, cf_status
import pytest

vars_to_test = ["tas","pr"]

input_convertors = [
    ("ERA5", {"variables": vars_to_test, "period": 2006, "freq": "hourly", "region": "europe"})
    # TODO: Add more datasets to test
]

@pytest.mark.parametrize("input_convertor", input_convertors)
def test_input_convertors(input_convertor):
    """
    Testing the input convertors through the input manager for the different datasets
    """
    manager = vp.InputManager(machine='hortense')
    ds = manager.load_data(input_convertor[0], **input_convertor[1])

    assert is_cf_compliant(ds), f"{cf_status(ds)}"