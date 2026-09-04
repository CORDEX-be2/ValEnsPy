from .diagnostic import (
    Diagnostic, Model2Self, Model2Ref, Ensemble2Ref, Ensemble2Self,
    match_ref_to_data, DEFAULT_IDENTITY_ATTRS,
)

# =============================================================================
# Pre-made diagnostics
# =============================================================================

from ._model2self import *
from ._model2ref import *
from ._ensemble2self import *
from ._ensemble2ref import *