import os
import sys
from pathlib import Path

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

import valenspy.input.converter_functions
from valenspy.input import InputConverter, INPUT_CONVERTORS
import valenspy.processing
from valenspy.diagnostic import Diagnostic, Model2Ref, Ensemble2Ref, Ensemble2Self
from valenspy.input import InputManager

from valenspy._utilities import is_cf_compliant, cf_status

# =============================================================================
# Demo-datafile
# =============================================================================

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

# demo files

demo_data_CF = os.path.join(
    BASE_PATH,
    "valenspy",
    "datafiles",
    "tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc",
)

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.1.0"
