from valenspy.diagnostic.diagnostic import Ensemble2Self
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
    "Ensemble_Quantile_Spatial_Mean",
    "Ensemble_Quantile_Closest_Member_Spatial_Mean"
]

Ensemble_Quantile_Spatial_Mean = Ensemble2Self(
    ensemble_quantile_of_spatial_mean,
    plot_map_per_dimension,
    "Ensemble quantiles of spatial mean",
    "The quantiles accross the ensembles spatial mean."
)

Ensemble_Quantile_Closest_Member_Spatial_Mean = Ensemble2Self(
    ensemble_quantile_closest_member_of_spatial_mean,
    plot_map,
    "Ensemble quantiles of closest member spatial mean",
    "The ensemble members that are closest to the quantiles of the spatial mean, and their spatial mean.",
    plot_type="facetted"
)