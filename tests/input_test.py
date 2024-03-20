import sys, os
from pathlib import Path

lib_folder = Path(__file__).resolve().parents[1].joinpath("src")

sys.path.insert(0, str(lib_folder))
import valenspy
import shapely

test_data_folder = Path(__file__).resolve().parents[0].joinpath("data")
file1 = test_data_folder.joinpath("tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc")
file2 = test_data_folder.joinpath("tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc")

#Testing the modeldata class

modeldata = valenspy.Modeldata(file_location=file1)
assert (modeldata.domain_bound == shapely.geometry.box(0, -89.4628215685774, 359.296875, 89.4628215685774)), "Domain is not correct"
assert (modeldata._is_CF_convention() == True), "Not in CF convention"

#Testing the ens_member class

ensmember = valenspy.Ensmember(data=[modeldata], experiment="historical", institution="EC-Earth3-Veg", model="EC-Earth3-Veg")

#Testing the inputprocessor class

## A dummy converter that does nothing - 
## The netCDF file is automatically converted to a Modeldata object and therefore CF convention is checked
simple_input = valenspy.InputProcessor(converter=lambda x: x)
assert (simple_input.convert_input(file1).file_location == file1), "Not returning the correct file"
assert ([md.file_location for md in simple_input.convert_input([file1, file2])] == [file1, file2]), "Not returning the correct files"

