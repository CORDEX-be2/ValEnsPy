import os
import sys
from pathlib import Path
from valenspy.input import InputConverter, INPUT_CONVERTORS
from valenspy.input import InputManager
#Processing
import valenspy.processing
#Diagnostic
from valenspy.diagnostic import Diagnostic, Model2Ref, Ensemble2Ref, Ensemble2Self
from valenspy.diagnostic.visualizations import (
    plot_diurnal_cycle,
    plot_time_series,
    plot_map,
    plot_spatial_bias,
    plot_maps_mod_ref_diff,
    plot_time_series_mod_ref,
    plot_points_on_map,
)
#Utility
from valenspy._utilities import is_cf_compliant, cf_status

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.1.0"
