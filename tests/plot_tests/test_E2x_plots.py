import pytest
import os
import xarray as xr
import matplotlib.pyplot as plt
import DataTree 

import valenspy

test_path = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(test_path, "plots")

ds = xr.tutorial.open_dataset('air_temperature').isel(time=slice(0, 10))

#Create an ensemble of 10 members.
datasets = {}
for i in range(10):
    datasets[f"member_{i}"] = ds + round(i/10, 2)

dt = DataTree.from_dict(datasets)

ds_mon = xr.tutorial.open_dataset("ersstv5")

ds_europe = ds.copy()
ds_europe['lat'] = ds_europe['lat'] + 10
ds_europe['lon'] = ds_europe['lon'] - 220


##########################
# Ensemble2Self Diagnostics #
##########################

m2s_diagnostics_tests = [
    (valenspy.diagnostic.DiurnalCycle, ds, "air"), 
    (valenspy.diagnostic.AnnualCycle, ds_mon, "sst"), 
    (valenspy.diagnostic.TimeSeriesSpatialMean, ds, "air")
    ]

@pytest.mark.parametrize("diagnostic,ds,var", m2s_diagnostics_tests)
def test_basic_M2S_plot(diagnostic, ds, var):
    result = diagnostic.apply(ds)
    diagnostic.plot(result[var])
    plt.savefig(os.path.join(output_path, f"M2S_{diagnostic.name}_basic.png"))
    plt.close()

@pytest.mark.parametrize("diagnostic,ds,var", m2s_diagnostics_tests)
def test_personalized_M2S_plot(diagnostic, ds, var):
    result = diagnostic.apply(ds)
    fig,ax = plt.subplots(figsize=(10, 5))
    diagnostic.plot(result[var], ax=ax, color="red", alpha=0.5)
    ax.set_xlabel("Personalized X-axis")
    ax.set_ylabel("Personalized Y-axis")
    plt.savefig(os.path.join(output_path, f"M2S_{diagnostic.name}_personalized.png"))
    plt.close()

m2s_europe_diagnostics_tests = [
    (valenspy.diagnostic.DiurnalCycle, ds_europe, "air"), 
    (valenspy.diagnostic.AnnualCycle, ds_mon, "sst"), 
    (valenspy.diagnostic.TimeSeriesSpatialMean, ds_europe, "air")
    ]

@pytest.mark.parametrize("diagnostic,ds,var", m2s_europe_diagnostics_tests)
def test_facetted_M2S_plot(diagnostic, ds, var):
    result = diagnostic.apply(ds, mask="prudence")
    diagnostic.plot(result[var], col="region", col_wrap=3)
    plt.savefig(os.path.join(output_path, f"M2S_{diagnostic.name}_facetted.png"))
    plt.close()

##########################
# Model2Ref Diagnostics ##
##########################

m2r_diagnostics_tests = [
    (valenspy.diagnostic.SpatialBias, ds, "air"),
    (valenspy.diagnostic.TemporalBias, ds, "air"),
    (valenspy.diagnostic.DiurnalCycleBias, ds, "air")
    ]

@pytest.mark.parametrize("diagnostic,ds,var", m2r_diagnostics_tests)
def test_basic_M2R_plot(diagnostic, ds, var):
    ref = ds + 2
    result = diagnostic.apply(ds, ref)
    diagnostic.plot(result[var])
    plt.savefig(os.path.join(output_path, f"M2R_{diagnostic.name}_basic.png"))
    plt.close()

@pytest.mark.parametrize("diagnostic,ds,var", m2r_diagnostics_tests)
def test_personalized_M2R_plot(diagnostic, ds, var):
    ref = ds + 2
    result = diagnostic.apply(ds, ref)
    fig,ax = plt.subplots(figsize=(10, 5))
    if len(result[var].dims) == 1:
        diagnostic.plot(result[var], ax=ax, color="red", alpha=0.5)
    else:
        diagnostic.plot(result[var], ax=ax, cmap="Blues", alpha=0.5, cbar_kwargs={'orientation': 'horizontal'})
    ax.set_xlabel("Personalized X-axis")
    ax.set_ylabel("Personalized Y-axis")
    plt.savefig(os.path.join(output_path, f"M2R_{diagnostic.name}_personalized.png"))
    plt.close()

m2r_diagnostics_tests_europe = [
    (valenspy.diagnostic.SpatialBias, ds_europe, "air"),
    (valenspy.diagnostic.TemporalBias, ds_europe, "air"),
    (valenspy.diagnostic.DiurnalCycleBias, ds_europe, "air")
    ]

@pytest.mark.parametrize("diagnostic,ds,var", m2r_diagnostics_tests_europe)
def test_facetted_M2R_plot(diagnostic, ds, var):
    ref = ds + 2
    result = diagnostic.apply(ds, ref, mask="prudence")
    diagnostic.plot(result[var], col="region", col_wrap=3)
    plt.savefig(os.path.join(output_path, f"M2R_{diagnostic.name}_facetted.png"))
    plt.close()