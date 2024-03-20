import src.valenspy
import xarray as xr
import dask 

ds = xr.open_dataset(r"C:\Users\kvandela\Projects\ValEnsPy\tests\data\tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc")

modeldata = src.valenspy.Modeldata(file_location=r"C:\Users\kvandela\Projects\ValEnsPy\tests\data\tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc")
print(modeldata.get_domain())

ensmember = src.valenspy.Ensmember(data=[modeldata, modeldata])
print(ensmember)


#Creating a simple input processor
simple_input = src.valenspy.InputProcessor(converter=lambda x: x)
modeldata = simple_input.convert_input(r"C:\Users\kvandela\Projects\ValEnsPy\tests\data\tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc")

#Using it for multiple files
data = simple_input.convert_input([r"C:\Users\kvandela\Projects\ValEnsPy\tests\data\tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc", 
                                   r"C:\Users\kvandela\Projects\ValEnsPy\tests\data\tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_195301-195312.nc"])

ensmember = src.valenspy.Ensmember(data=data, experiment="historical", institution="EC-Earth3-Veg", model="EC-Earth3-Veg")
ensmember

