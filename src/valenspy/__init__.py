import os
import sys
from pathlib import Path
import valenspy.input.converter_functions
from valenspy.input import InputConverter, INPUT_CONVERTORS
import valenspy.processing
from valenspy.diagnostic import Diagnostic, Model2Ref, Ensemble2Ref, Ensemble2Self
from valenspy.input import InputManager

from valenspy._utilities import is_cf_compliant, cf_status

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.1.0"
