from valenspy.diagnostic.diagnostic import Ensemble2Self
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

import seaborn as sns

__all__ = [
    "Ensemble_Histogram",
    "Ensemble_Quantile_Spatial_Mean",
    "Ensemble_Quantile_Closest_Member_Spatial_Mean",
    "Ensemble_Quantile_Of_Climate_Change_Signal",
    "Ensemble_Quantile_Closest_Member_Of_Climate_Change_Signal",
    "Ensemble_Spatial_Mean",
    "ClimateChangeSignalPerMember",
    "ClimatologyPerMember",
    "ClimateChangeSignalEnsembleMeanGrid",
    ]
Ensemble_Histogram = Ensemble2Self(
    ensemble_member_means,
    sns.histplot,
    "Ensemble member means",
    "The histogram of the ensemble member means."
)

Ensemble_Quantile_Spatial_Mean = Ensemble2Self(
    ensemble_quantile_of_spatial_mean,
    plot_quantile_map,
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

Ensemble_Quantile_Of_Climate_Change_Signal = Ensemble2Self(
    ensemble_quantile_of_climate_change_signal,
    plot_quantile_map,
    "Ensemble quantiles of climate change signal",
    "The quantiles across the ensemble's climate-change signal (future period minus "
    "historical, matched per member by identity) - the per-member spread of the "
    "projected change at one future period, summarized by quantile. Call as "
    "Ensemble_Quantile_Of_Climate_Change_Signal(dt, historical=\"historical\", "
    "period=\"ssp245\", quantile=[0.1, 0.5, 0.9])."
)

Ensemble_Quantile_Closest_Member_Of_Climate_Change_Signal = Ensemble2Self(
    ensemble_quantile_closest_member_of_climate_change_signal,
    plot_map,
    "Ensemble quantiles of closest member climate change signal",
    "The ensemble members whose climate-change signal is closest to each requested "
    "quantile, and their own signal - the closest-member counterpart to "
    "Ensemble_Quantile_Of_Climate_Change_Signal. Call as "
    "Ensemble_Quantile_Closest_Member_Of_Climate_Change_Signal(dt, historical=\"historical\", "
    "period=\"ssp245\", quantile=[0.1, 0.5, 0.9], var=\"tas\").",
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
    "together. Call as ClimateChangeSignalPerMember(dt, historical=\"historical\", "
    "future_periods=[\"ssp245\", \"ssp585\"])."
)

ClimatologyPerMember = Ensemble2Self(
    climatology_per_member,
    plot_reference_future_periods_grid,
    "Climatology per member",
    "For each ensemble member individually, its reference-period mean followed by its own time "
    "mean for each future period (not a change signal - see ClimateChangeSignalPerMember for "
    "that) - unlike Ensemble_Spatial_Mean, members are kept separate rather than averaged "
    "together. Call as ClimatologyPerMember(dt, historical=\"historical\", "
    "future_periods=[\"ssp245\", \"ssp585\"])."
)

ClimateChangeSignalEnsembleMeanGrid = Ensemble2Self(
    climate_change_signal_ensemble_mean_grid,
    plot_reference_future_periods_grid,
    "Climate change signal of the ensemble mean",
    "The reference-period ensemble mean next to the ensemble-mean climate change signal for "
    "each future period - the ensemble-mean counterpart to ClimateChangeSignalPerMember "
    "(members averaged together here, rather than kept as separate rows). Call as "
    "ClimateChangeSignalEnsembleMeanGrid(dt, historical=\"historical\", "
    "future_periods=[\"ssp245\", \"ssp585\"])."
)
