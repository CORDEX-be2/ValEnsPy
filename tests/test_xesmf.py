import xesmf as xe
import xarray as xr
import numpy as np
import valenspy as vp

# Load the data
ds = xr.open_dataset("/dodrio/scratch/projects/2022_200/external/era5/europe/2m_temperature/monthly/2m_temperature_era5_europe_monthly_min_2023.nc", chunks="auto")
ds

# Load the target grid
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(16, 75, 1.0), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(200, 330, 1.5), {"units": "degrees_east"}),
    }
)
ds_out

# # Create the regridder
# regridder = xe.Regridder(ds, ds_out, "bilinear")
# ds_reg = regridder(ds)

ds_reg = vp.processing.regrid.remap_xesmf(ds, ds_out, method="bilinear", regridding_kwargs={"keep_attrs":True})