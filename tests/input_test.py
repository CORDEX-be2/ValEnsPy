import sys, os
from pathlib import Path
import numpy

lib_folder = Path(__file__).resolve().parents[1].joinpath("src")

sys.path.insert(0, str(lib_folder))
import valenspy
import shapely

test_data_folder = Path(__file__).resolve().parents[0].joinpath("data")
file1 = test_data_folder.joinpath("tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc")
file2 = test_data_folder.joinpath("tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195401-195412.nc")

#Testing the modeldata class

md = valenspy.Modeldata(file_location=file1)
assert (md.domain_bound == shapely.geometry.box(0, -89.4628215685774, 359.296875, 89.4628215685774)), "Domain is not correct"
assert (md._is_CF_convention() == True), "Not in CF convention"

#Testing the ens_member class

mem = valenspy.Ensmember(data=[md, valenspy.Modeldata(file_location=file2)], experiment="historical", institution="EC-Earth3-Veg", model="EC-Earth3-Veg")
assert (mem.time_period == (numpy.datetime64('1953-01-16T12:00:00.000000000'), numpy.datetime64('1954-12-16T12:00:00.000000000'))), "Time period is not correct"
assert (mem.get_domain == shapely.geometry.box(0, -89.4628215685774, 359.296875, 89.4628215685774)), "Domain is not correct"
assert (mem.resolution == "100 km"), "Resolution is not correct"

#Testing the ensemble class

#Testing the inputprocessor class

## A dummy converter that does nothing - 
## The netCDF file is automatically converted to a Modeldata object and therefore CF convention is checked
simple_input = valenspy.InputProcessor(converter=lambda x: x)
assert (simple_input.convert_input(file1).file_location == file1), "Not returning the correct file"
assert ([md.file_location for md in simple_input.convert_input([file1, file2])] == [file1, file2]), "Not returning the correct files"
