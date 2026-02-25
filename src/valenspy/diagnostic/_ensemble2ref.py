from valenspy.diagnostic.diagnostic import Ensemble2Ref
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
    "ClimateChangeSignalOfSpatialMean",
    "MetricsRankings"
    ]

ClimateChangeSignalOfSpatialMean = Ensemble2Ref(
    climate_change_signal_of_spatial_mean,
    plot_map,
    "Climate Change Signal of the spatial means",
    "The spatial climate change signal as the difference between the temporal average of two periods",
    plot_type="facetted"
)

ClimateChangeSignal = Ensemble2Ref(
    mean_climate_change_signal,
    None,
    "Climate Change Signal",
    "The climate change signal as the difference between the spatial and temporal average of two periods."
)

MetricsRankings = Ensemble2Ref(
    calc_metrics_dt,
    plot_metric_ranking,
    "Metrics Rankings",
    "The rankings of ensemble members with respect to several metrics when compared to the reference."
)