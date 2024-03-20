import os
import sys
from pathlib import Path

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

from valenspy.modeldata import Modeldata
from valenspy.ensmember import Ensmember
from valenspy.ensemble import Ensemble
from valenspy.inputprocessor import InputProcessor

    
# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.0.0"