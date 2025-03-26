import xclim
from xclim.testing import open_dataset
import xarray as xr
from datatree import DataTree
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
        
