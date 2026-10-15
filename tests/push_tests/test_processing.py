import xclim
from xclim.testing import open_dataset
import xarray as xr
from xarray import DataTree
import valenspy as vp

def test_unit_conversion():
    #xclim test dataset
    ds = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")
    dt = DataTree.from_dict({"data": ds})

    conversions = {
        "tasmax": "°C",
        "pr": "mm d-1"
    }

    ds_c = ds.copy()
    dt_c = dt.copy()
    for var, target_unit in conversions.items():
        ds_c = vp.convert_units_to(ds_c, var, target_unit)
        dt_c = vp.convert_units_to(dt_c, var, target_unit)

        assert ds_c[var].units == target_unit, f"Unit conversion failed for xr.Dataset variable {var}: was {ds_c[var].units} expected {target_unit}"
        assert dt_c.data[var].units == target_unit, f"Unit conversion failed for DataTree variable {var}: was {dt_c.data[var].units} expected {target_unit}"

    assert isinstance(ds_c, xr.Dataset), "xr.Dataset not returned after unit conversion"
    assert isinstance(dt_c, DataTree), "DataTree not returned after unit conversion"
        
def test_convert_units_by_mapping():
    ds = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")
    dt = DataTree.from_dict({"a": ds, "b": ds})

    units_map = {"K": "degC", "kg m-2 s-1": ("mm/d", "hydro")}
    ds_c = vp.convert_units_by_mapping(ds, units_map)
    dt_c = vp.convert_units_by_mapping(dt, units_map)

    for out in (ds_c, dt_c["a"].to_dataset(), dt_c["b"].to_dataset()):
        assert out["tasmax"].units == "°C"
        assert out["pr"].units == "mm d-1"
    assert isinstance(ds_c, xr.Dataset) and isinstance(dt_c, DataTree)
    assert ds["tasmax"].units == "K", "input must not be modified"

    custom = vp.convert_units_by_mapping(ds, {"K": "degF"})
    assert custom["tasmax"].units == "°F" and custom["pr"].units == ds["pr"].units

def test_convert_units_by_mapping_temperature_difference():
    diff = xr.DataArray([6.0], dims="x", attrs={"units": "K", "units_metadata": "temperature: difference"})
    absolute = xr.DataArray([300.0], dims="x", attrs={"units": "K"})
    ds = xr.Dataset({"excess": diff, "tas": absolute})
    dt = DataTree.from_dict({"a": ds, "b": xr.Dataset({"tas": absolute})})

    for out in (vp.convert_units_by_mapping(ds, {"K": "degC"}),
                vp.convert_units_by_mapping(dt, {"K": "degC"})["a"].to_dataset()):
        assert float(out["excess"]) == 6.0, "a difference must not get the 273.15 offset"
        assert out["excess"].attrs["units_metadata"] == "temperature: difference"
        assert abs(float(out["tas"]) - 26.85) < 1e-6

    assert abs(float(vp.convert_units_by_mapping(ds, {"K": "degF"})["excess"]) - 10.8) < 1e-6

def test_convert_units_by_mapping_label():
    dd = xr.DataArray([100.0], dims="x", attrs={"units": "K days"})
    ds = xr.Dataset({"dd": dd, "other": dd.copy()})
    dt = DataTree.from_dict({"a": ds, "b": xr.Dataset({"dd": dd})})
    units_map = {"K days": {"units": "degC days", "label": "°C days"}}

    out = vp.convert_units_by_mapping(ds, units_map)
    assert float(out["dd"]) == 100.0 and out["dd"].units == "°C days"
    out_dt = vp.convert_units_by_mapping(dt, units_map)
    assert out_dt["a"].to_dataset()["dd"].units == "°C days" and out_dt["b"].to_dataset()["dd"].units == "°C days"
    assert vp.convert_units_by_mapping(ds, {"K days": {"units": "degC days"}})["dd"].units == "d degC"

def test_xclim_indicators():
    #xclim test dataset
    ds = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")
    dt = DataTree.from_dict({"data": ds, "data2": (ds+ds)/1.9})

    #xclim indicator
    indicators = {
        "bedd" : {"indicator": xclim.indicators.atmos.biologically_effective_degree_days, "vars": ["tasmin", "tasmax", "lat"], "units": "K days"},
        "txx" : {"indicator": xclim.indicators.cf.txx, "vars":"tasmax", "units": "degree_Celsius"},
    }

    for name, ind in indicators.items():
        dt_i = vp.xclim_indicator(dt, ind["indicator"], ind["vars"])

        assert name in dt_i.data, f"Indicator {name} not found in the output DataTree"
        assert ind["units"] in dt_i.data2[name].units, f"Units {ind['units']} not found in the output DataTree"

    assert isinstance(dt_i, DataTree), "DataTree not returned after indicator calculation"
