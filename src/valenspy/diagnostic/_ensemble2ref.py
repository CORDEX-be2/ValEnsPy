from valenspy.diagnostic.diagnostic import Ensemble2Ref
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
    "ClimateChangeSignalOfSpatialMean",
    "ClimateChangeSignalEnsembleMean",
    "ClimateChangeSignal",
    "MetricsRankings",
    ]

ClimateChangeSignalOfSpatialMean = Ensemble2Ref(
    climate_change_signal_of_spatial_mean,
    plot_map,
    "Climate Change Signal of the spatial means",
    "The spatial climate change signal as the difference between the temporal average of two periods",
    plot_type="facetted"
)

ClimateChangeSignalEnsembleMean = Ensemble2Ref(
    climate_change_signal_ensemble_mean,
    plot_ensemble_mean_map,
    "Climate Change Signal of the ensemble means",
    "The ensemble climate change signal as the difference between the temporal average of two periods",
)

ClimateChangeSignal = Ensemble2Ref(
    mean_climate_change_signal,
    lambda ds : ds,
    "Climate Change Signal",
    "The climate change signal as the difference between the spatial and temporal average of two periods."
)

MetricsRankings = Ensemble2Ref(
    calc_metrics_dt,
    plot_metric_ranking,
    "Metrics Rankings",
    "The rankings of ensemble members with respect to several metrics when compared to the reference."
)