import pytest

def test_access_main():
    """
    Test to see if the main functions are accessible from the main package.
    """
    import valenspy as vp

    #Input related functions/Objects - directly accessible
    assert vp.INPUT_CONVERTORS, "INPUT_CONVERTORS not accessible from valenspy"
    assert vp.InputConverter, "InputConverter not accessible from valenspy"
    assert vp.InputManager, "InputManager not accessible from valenspy"

    #Preprocessing - directly accessible

    #Diagnostic related functions/Objects - directly accessible
    ##Objects
    assert vp.Diagnostic, "Diagnostic not accessible from valenspy"
    assert vp.Model2Ref, "Model2Ref not accessible from valenspy"
    assert vp.Ensemble2Ref, "Ensemble2Ref not accessible from valenspy"
    assert vp.Ensemble2Self, "Ensemble2Self not accessible from valenspy"

    ##Visualizations functions
    assert vp.plot_diurnal_cycle, "plot_diurnal_cycle not accessible from valenspy"
    assert vp.plot_time_series, "plot_time_series not accessible from valenspy"
    assert vp.plot_map, "plot_map not accessible from valenspy"
    assert vp.plot_spatial_bias, "plot_spatial_bias not accessible from valenspy"
    assert vp.plot_maps_mod_ref_diff, "plot_maps_mod_ref_diff not accessible from valenspy"
    assert vp.plot_time_series_mod_ref, "plot_time_series_mod_ref not accessible from valenspy"
    assert vp.plot_points_on_map, "plot_points_on_map not accessible from valenspy"

def test_access_input():
    """
    Test to see if the input related functions are accessible from the valenspy.input module.
    """
    import valenspy.input  

    #Input converter functions
    assert valenspy.input.EOBS_to_CF, "EOBS_to_CF converter function not accessible from valenspy.input"
    assert valenspy.input.ERA5_to_CF, "ERA5_to_CF converter function not accessible from valenspy.input"
    assert valenspy.input.ERA5Land_to_CF, "ERA5Land_to_CF converter function not accessible from valenspy.input"
    assert valenspy.input.CLIMATE_GRID_to_CF, "CLIMATE_GRID_to_CF converter function not accessible from valenspy.input"
    assert valenspy.input.CCLM_to_CF, "CCLM_to_CF converter function not accessible from valenspy.input"
    assert valenspy.input.ALARO_K_to_CF, "ALARO_K_to_CF converter function not accessible from valenspy.input"
    assert valenspy.input.RADCLIM_to_CF, "RADCLIM_to_CF converter function not accessible from valenspy.input"

def test_access_diagnostics():
    """
    Test to see if the diagnostic related functions are accessible from the valenspy.diagnostic module.
    """
    import valenspy.diagnostic

    #Pre-made diagnostics
    assert valenspy.diagnostic.DiurnalCycle, "DiurnalCycle not accessible from valenspy.diagnostic"
    assert valenspy.diagnostic.TimeSeriesSpatialMean, "TimeSeriesSpatialMean not accessible from valenspy.diagnostic"
    assert valenspy.diagnostic.SpatialBias, "SpatialBias not accessible from valenspy.diagnostic"
    assert valenspy.diagnostic.TemporalBias, "TemporalBias not accessible from valenspy.diagnostic"
    assert valenspy.diagnostic.DiurnalCycle, "DiurnalCycle not accessible from valenspy.diagnostic"