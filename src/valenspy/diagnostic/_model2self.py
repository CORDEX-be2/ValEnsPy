from valenspy.diagnostic.diagnostic import Model2Self
from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

__all__ = [
    "DiurnalCycle",
    "AnnualCycle",
    "TimeSeriesSpatialMean",
    "TimeSeriesTrendSpatialMean",
    "SpatialTimeMean",
    "UrbanHeatIsland",
    "UrbanHeatIslandDiurnalCycle"
]

DiurnalCycle = Model2Self(
    diurnal_cycle, 
    plot_diurnal_cycle, 
    "Diurnal Cycle", 
    "The diurnal cycle of the data.",
    plot_type="single"
)
AnnualCycle = Model2Self(
    annual_cycle,
    plot_annual_cycle,
    "Annual Cycle",
    "The annual cycle of the data.",
    plot_type="single"
)
TimeSeriesSpatialMean = Model2Self(
    time_series_spatial_mean,
    plot_time_series,
    "Time Series",
    "The time series of the data - if the data is spatial, the spatial mean is taken.",
    plot_type="single"
)
TimeSeriesTrendSpatialMean = Model2Self(
    time_series_trend,
    plot_time_series,
    "Time Series Trend",
    "The time series trend of the data - if the data is spatial, the spatial mean is taken.",
    plot_type="single"
)
SpatialTimeMean = Model2Self(
    spatial_time_mean,
    plot_map,
    "Spatial Mean",
    "The spatial representation of the time mean of the data."
)
UrbanHeatIsland = Model2Self(
    urban_heat_island,
    plot_time_series,
    "Urban Heat Island",
    "The urban heat island as the difference in temperature between urban and rural areas.",
)
UrbanHeatIslandDiurnalCycle = Model2Self(
    urban_heat_island_diurnal_cycle,
    plot_diurnal_cycle,
    "Urban Heat Island Diurnal Cycle",
    "The diurnal cycle of the urban heat island.",
)