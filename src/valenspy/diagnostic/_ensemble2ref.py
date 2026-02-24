from valenspy.diagnostic.diagnostic import Ensemble2Ref
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
    "ClimateChangeSignal",
    "MetricsRankings"
    ]

ClimateChangeSignal = Ensemble2Ref(
    climate_change_signal,
    plot_map,
    "Climate Change Signal",
    "The climate change signal as the average difference between the GWL and the reference period.",
    plot_type="facetted"
)

MetricsRankings = Ensemble2Ref(
    calc_metrics_dt,
    plot_metric_ranking,
    "Metrics Rankings",
    "The rankings of ensemble members with respect to several metrics when compared to the reference."
)