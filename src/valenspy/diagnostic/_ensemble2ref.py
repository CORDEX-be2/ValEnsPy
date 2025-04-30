from valenspy.diagnostic.diagnostic import Ensemble2Ref
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = ["MetricsRankings"]

MetricsRankings = Ensemble2Ref(
    calc_metrics_dt,
    plot_metric_ranking,
    "Metrics Rankings",
    "The rankings of ensemble members with respect to several metrics when compared to the reference."
)