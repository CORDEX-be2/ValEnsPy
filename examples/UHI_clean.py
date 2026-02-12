# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
from dask.diagnostics import ProgressBar
from pathlib import Path
import xarray as xr
import matplotlib.pyplot as plt
from datatree import DataTree

import valenspy as vp

# %%
ds = xr.open_mfdataset("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/CORDEXbeII/run_ALARO_sfx/out/PROD_ERA5_NEST_BE4_old/19780101/CORDEX/CMIP6/DD/BE-04/RMIB-UGent/ERA5-ALARO1-SFX/evaluation/r1i1p1f1/ALARO1-SFX/v1-r1/1hr/tas/v20250121/*.nc", chunks="auto", decode_coords="all")
ds2 = xr.open_mfdataset("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/CORDEXbeII/run_ALARO_sfx/out/PROD_ERA5_NEST_BE4_1/19780101/CORDEX/CMIP6/DD/BE-04/RMIB-UGent/ERA5-ALARO1-SFX/evaluation/r1i1p1f1/ALARO1-SFX/v1-r1/1hr/tas/v20250121/*.nc", chunks="auto", decode_coords="all")
ds3 = xr.open_mfdataset("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/CORDEXbeII/run_ALARO_sfx/out/PROD_ERA5_NEST_BE4_1_dry/19780101/CORDEX/CMIP6/DD/BE-04/RMIB-UGent/ERA5-ALARO1-SFX/evaluation/r1i1p1f1/ALARO1-SFX/v1-r1/1hr/tas/v20250121/*.nc", chunks="auto", decode_coords="all")

# %%
dt = DataTree.from_dict({
    "control": ds,
    "mlch_hum" : ds2,
    "mlch_dry" : ds3
})

# %%
from valenspy.diagnostic import UrbanHeatIslandDiurnalCycle

melle = (50.980438, 3.815763) #lat, lon
gent_center = (51.04573, 3.72449) #lat, lon

# %%
dt_sel = dt.sel(time=slice("1978-01-01", "1978-01-31"))
with ProgressBar():
    dt_uhi = UrbanHeatIslandDiurnalCycle(dt_sel, urban_coord=gent_center, rural_coord=melle, projection="lcc")
    dt_uhi = dt_uhi.compute()

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
UrbanHeatIslandDiurnalCycle.plot_dt(dt_uhi, ax=ax, var="tas")
ax.legend()
plt.title("Urban Heat Island Diurnal Cycle - January 1978")
plt.show()
