import os
import sys
from pathlib import Path

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

from valenspy.inputconverter import InputConverter
from valenspy.preprocessor import Preprocessor
from valenspy.diagnostic import Diagnostic, Model2Ref, Ensemble2Ref
import valenspy.preprocessing_tasks
import valenspy.cf_checks

# =============================================================================
# Demo-datafile
# =============================================================================

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

# demo files

demo_data_CF = os.path.join(
    BASE_PATH, "valenspy", "datafiles", "tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc"
)

# =============================================================================
# Pre-made diagnostics
# =============================================================================

# Load some pre-made diagnostics which others can use and contribute to

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.1.0"
