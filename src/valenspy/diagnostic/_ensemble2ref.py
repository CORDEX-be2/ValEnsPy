from valenspy.diagnostic.diagnostic import Ensemble2Ref
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
        "MetricsRankings",
        "EnsembleSubSelection",
        ]

MetricsRankings = Ensemble2Ref(
    calc_metrics_dt,
    plot_metric_ranking,
    "Metrics Rankings",
    "The rankings of ensemble members with respect to several metrics when compared to the reference."
)

EnsembleSubSelection = Ensemble2Ref(
    case_sub_selection,
    {"default":
        default_plot_kwargs({
        "x": "var",
        "y": "abs_change",
        "selected": ["highest", "middle", "lowest"],
        "sel_colors": {"highest": "red", "middle": "blue", "lowest": "green"}
            })(ensemble_selection_boxplot),
    "heatmap":
        default_plot_kwargs({
        "index": "label",
        "columns": "var",
        "values": "rel_change"
        })(ensemble_change_signal_heatmap)},
    "Ensemble Sub Selection",
    "The sub selection of ensemble members."
    )
