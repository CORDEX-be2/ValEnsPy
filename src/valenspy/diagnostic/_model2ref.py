from valenspy.diagnostic.diagnostic import Model2Ref
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
    "SpatialBias",
    "TemporalBias",
    "DiurnalCycleBias"
]

# Model2Ref diagnostics
SpatialBias = Model2Ref(
    spatial_bias,
    plot_map,
    "Spatial Bias",
    "The spatial bias of the data compared to the reference.",
    plot_type="facetted"
)
TemporalBias = Model2Ref(
    temporal_bias,
    plot_time_series,
    "Temporal Bias",
    "The temporal bias of the data compared to the reference.",
    plot_type="single"
)
DiurnalCycleBias = Model2Ref(
    diurnal_cycle_bias,
    plot_diurnal_cycle,
    "Diurnal Cycle Bias",
    "The diurnal cycle bias of the data compared to the reference.",
    plot_type="single"
)