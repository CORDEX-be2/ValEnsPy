from valenspy.diagnostic.diagnostic import Ensemble2Self
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

import seaborn as sns

__all__ = [
    "Ensemble_Histogram",
    "Ensemble_Quantile_Spatial_Mean",
    "Ensemble_Quantile_Closest_Member_Spatial_Mean",
    "Ensemble_Spatial_Mean",
    "ClimateChangeSignalPerMember",
    "ClimatologyPerMember",
    ]
Ensemble_Histogram = Ensemble2Self(
    ensemble_member_means,
    sns.histplot,
    "Ensemble member means",
    "The histogram of the ensemble member means."
)

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

Ensemble_Spatial_Mean = Ensemble2Self(
    ensemble_spatial_mean,
    plot_map,
    "Ensemble spatial mean",
    "The spatial mean across the ensemble members.",
)

ClimateChangeSignalPerMember = Ensemble2Self(
    climate_change_signal_per_member,
    plot_reference_future_periods_grid,
    "Climate change signal per member",
    "For each ensemble member individually, its reference-period mean followed by its own "
    "climate change signal (future minus its own reference period) for each future period - "
    "unlike ClimateChangeSignalEnsembleMean, members are kept separate rather than averaged "
    "together. Call as ClimateChangeSignalPerMember(ref, fut_periods={{label: DataTree, ...}})."
)

ClimatologyPerMember = Ensemble2Self(
    climatology_per_member,
    plot_reference_future_periods_grid,
    "Climatology per member",
    "For each ensemble member individually, its reference-period mean followed by its own time "
    "mean for each future period (not a change signal - see ClimateChangeSignalPerMember for "
    "that) - unlike Ensemble_Spatial_Mean, members are kept separate rather than averaged "
    "together. Call as ClimatologyPerMember(ref, fut_periods={{label: DataTree, ...}})."
)
