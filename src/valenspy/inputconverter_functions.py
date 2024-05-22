from pathlib import Path
import xarray as xr
from yaml import safe_load

# get current path
files = Path(__file__).resolve().parent 

# open CORDEX variable lookup dictionary
with open(files / "ancilliary_data" / "CORDEX_variables.yml") as file:
    CORDEX_VARIABLES = safe_load(file)



def _non_convertor(file: Path) -> Path:
    """A dummy function that does not convert the input file as it already is in the correct format."""
    return file



def EOBS_to_CF(file: Path) -> Path:

    """
    Convert the EOBS netCDF file to a netCDF file in CF convention 

    Parameters
    ----------
    file : Path
        The path to the netCDF file of specific variable to convert 
    
    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """

    # open the observation dataset
    ds = xr.open_mfdataset(files, combine='by_coords', chunks='auto')

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator? 
    obsdata_name = "EOBS"

    with open(files / "ancilliary_data" / Path(obsdata_name+"_lookup.yml")) as file:
        obs_LOOKUP = safe_load(file)

    
    # make EOBS CF compliant

    # update variable name 
    ds = ds.rename_vars({obs_var: variable})
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})

    # Unit conversion - hard coded EOBS units for units different to CORDEX
    if obs_LOOKUP[variable]['obs_units'] == 'Celcius': 
        ds[variable] = ds[variable] + 273.15 # Celcius to Kelvin
    elif obs_LOOPUP[variable]['obs_units'] == 'hPa': 
        ds[variable] = ds[variable] * 100 # hPa to Pa
    elif obs_LOOPUP[variable]['obs_units'] == 'mm': # ! note observations remain daily time frequency
        ds[variable] = ds[variable] / 86400 # mm to kg m^-2 s^-1 

    # update unit attribute
    ds.attrs["units"] = CORDEX_VARIABLES[variable]["units"] # from the CORDEX look-up table 

    # add necessary metadata
    ds.attrs["standard_name"]      = CORDEX_VARIABLES[variable]["standard_name"]  # from the CORDEX look-up table
    ds.attrs["long_name"]          = CORDEX_VARIABLES[variable]["long_name"]  # from the CORDEX look-up table
    ds.attrs["original_name"]      = obs_LOOKUP[variable]["obs_name"]
    ds.attrs["original_long_name"] = obs_LOOKUP[variable]["obs_long_name"]

    # additional attributes -- hard coded for EOBS
    ds.attrs["time_freq"]          = "daily"  # possible values: daily, hourly, monthly, yearly
    ds.attrs["spatial_resolution"] = "0.1deg"  
    ds.attrs["domain"]             = "europe"

    ## question: check for CF compliance to be done here? 

    return ds
