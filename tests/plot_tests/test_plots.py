import pytest
import os
import xarray as xr
import matplotlib.pyplot as plt

import valenspy

test_path = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(test_path, "plots")

ds = xr.tutorial.open_dataset('air_temperature').isel(time=slice(0, 10))

##########################
# Model2Self Diagnostics #
##########################

m2s_diagnostics = [valenspy.diagnostic.DiurnalCycle, valenspy.diagnostic.TimeSeriesSpatialMean]

@pytest.mark.parametrize("diagnostic", m2s_diagnostics)
@pytest.mark.parametrize("ds", [ds])
def test_basic_M2S_plot(diagnostic, ds):
    result = diagnostic.apply(ds)
    diagnostic.plot(result.air)
    plt.savefig(os.path.join(output_path, f"M2S_{diagnostic.name}_basic.png"))
    plt.close()

@pytest.mark.parametrize("diagnostic", m2s_diagnostics)
@pytest.mark.parametrize("ds", [ds])
def test_personalized_M2S_plot(diagnostic, ds):
    result = diagnostic.apply(ds)
    fig,ax = plt.subplots(figsize=(10, 5))
    diagnostic.plot(result.air, ax=ax, color="red", alpha=0.5)
    ax.set_xlabel("Personalized X-axis")
    ax.set_ylabel("Personalized Y-axis")
    plt.savefig(os.path.join(output_path, f"M2S_{diagnostic.name}_personalized.png"))
    plt.close()

##########################
# Model2Ref Diagnostics ##
##########################

m2r_diagnostics = [valenspy.diagnostic.SpatialBias, valenspy.diagnostic.TemporalBias, valenspy.diagnostic.DiurnalCycleBias]

@pytest.mark.parametrize("diagnostic", m2r_diagnostics)
@pytest.mark.parametrize("ds", [ds])
def test_basic_M2R_plot(diagnostic, ds):
    ref = ds + 2
    result = diagnostic.apply(ds, ref)
    diagnostic.plot(result.air)
    plt.savefig(os.path.join(output_path, f"M2R_{diagnostic.name}_basic.png"))
    plt.close()

@pytest.mark.parametrize("diagnostic", m2r_diagnostics)
@pytest.mark.parametrize("ds", [ds])
def test_personalized_M2R_plot(diagnostic, ds):
    ref = ds + 2
    result = diagnostic.apply(ds, ref)
    fig,ax = plt.subplots(figsize=(10, 5))
    if len(result.air.dims) == 1:
        diagnostic.plot(result.air, ax=ax, color="red", alpha=0.5)
    else:
        diagnostic.plot(result.air, ax=ax, cmap="blues", alpha=0.5, cbar_kwargs={'orientation': 'horizontal'})
    ax.set_xlabel("Personalized X-axis")
    ax.set_ylabel("Personalized Y-axis")
    plt.savefig(os.path.join(output_path, f"M2R_{diagnostic.name}_personalized.png"))
    plt.close()