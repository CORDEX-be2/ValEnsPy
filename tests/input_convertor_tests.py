from pathlib import Path

import xarray as xr
import pandas as pd

# The Valenspy package
import valenspy as vp
from valenspy.inputconverter_functions import EOBS_to_CF
from valenspy._utilities import load_yml
from valenspy.cf_checks import is_cf_compliant

manager = vp.InputManager(machine='hortense')

vars_to_test = ["tas","pr"]    

# ERA5 data
ds_era5 = manager.load_data("ERA5",vars_to_test, period=2000,freq="daily", region="europe", path_identifiers=["min"])
assert is_cf_compliant(ds_era5)

# EOBS data
ds_eobs = manager.load_data("EOBS",vars_to_test, path_identifiers=["_mean"])
assert is_cf_compliant(ds_eobs)
