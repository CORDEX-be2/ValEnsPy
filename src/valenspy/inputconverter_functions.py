from pathlib import Path
import xarray as xr

def _non_convertor(file: Path) -> Path:
    """A dummy function that does not convert the input file as it already is in the correct format."""
    return file

def EOBS_to_CF(file: Path) -> Path:
    """Convert the EOBS netCDF file to a netCDF file in CF convention."""
    #Open the EOBS file
    ds = xr.open_dataset(file)
    
    #Rename the dimensions
    ds = ds.rename_dims({"latitude": "lat", "longitude": "lon"})
    




