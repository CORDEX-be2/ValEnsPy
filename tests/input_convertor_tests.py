from pathlib import Path

import xarray as xr
import pandas as pd

# The Valenspy package
import valenspy as vp
from valenspy.inputconverter_functions import EOBS_to_CF
from valenspy._utilities import load_yml

# EOBs data
