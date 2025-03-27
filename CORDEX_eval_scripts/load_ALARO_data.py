import xarray as xr
import valenspy as vp
from datatree import DataTree
from pathlib import Path
import intake

cat = intake.open_esm_datastore("https://github.com/CORDEX-be2/ValEnsPy/raw/refs/heads/CORDEX_basic_eval/CORDEX_eval_scripts/CORDEX-CMIP6.json")

cat_subset = cat.search(frequency="day", variable_id=['tas', 'pr'])
data_set_dict = cat_subset.to_dataset_dict(xarray_open_kwargs={"decode_coords":"all"})

#Save the entry to a netcdf file
for key in data_set_dict.keys():
    data_set_dict[key].to_netcdf(f"/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/notebooks/intermediate_data/{key}.nc")

ds = xr.open_dataset("/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/vsc46032_kobe/ValEnsPy/notebooks/intermediate_data/CMIP6.CORDEX.EUR-11.CNRM-CERFACS.CNRM-CM6-1.historical.r1i1p1f2.day.pr.v20190219.nc")

