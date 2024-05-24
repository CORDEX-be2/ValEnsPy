from pathlib import Path
import xarray as xr
from yaml import safe_load
from valenspy.cf_checks import is_cf_compliant, cf_status

# get path of source code (current path)
src_path = Path(__file__).resolve().parent 

# open CORDEX variable lookup dictionary
with open(src_path / "ancilliary_data" / "CORDEX_variables.yml") as file:
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
    ds = xr.open_mfdataset(file, combine='by_coords', chunks='auto')

    # open observational specific lookyp dictionary - now hardcoded for EOBS, but this can be automated, potentially in the Path generator? 
    obsdata_name = "EOBS"

    with open(src_path / "ancilliary_data" / Path(obsdata_name+"_lookup.yml")) as file:
        obs_LOOKUP = safe_load(file)

    
    # make EOBS CF compliant

    for var_obs in ds.data_vars: 

        # Get the CORDEX variable in the observational dataset using the observational lookup table  
        var = next((k for k, v in obs_LOOKUP.items() if v.get('obs_name') == var_obs), None)

        # update variable name to CORDEX variable name
        ds = ds.rename_vars({obs_LOOKUP[var]["obs_name"]: var})

        # from here on, use CORDEX variable name to access data array and do rest of conversion

        # Unit conversion - hard coded EOBS units for units different to CORDEX
        if obs_LOOKUP[var]['obs_units'] == 'Celcius': 
            ds[var] = ds[var] + 273.15 # Celcius to Kelvin
        elif obs_LOOKPUP[var]['obs_units'] == 'hPa': 
            ds[var] = ds[var] * 100 # hPa to Pa
        elif obs_LOOKPUP[var]['obs_units'] == 'mm': # ! note observations remain daily time frequency
            ds[var] = ds[var] / 86400 # mm to kg m^-2 s^-1 

        # update unit attribute
        ds[var].attrs["units"] = CORDEX_VARIABLES[var]["units"] # from the CORDEX look-up table 

        # add necessary metadata
        ds[var].attrs["standard_name"]      = CORDEX_VARIABLES[var]["standard_name"]  # from the CORDEX look-up table
        ds[var].attrs["long_name"]          = CORDEX_VARIABLES[var]["long_name"]  # from the CORDEX look-up table
        ds[var].attrs["original_name"]      = obs_LOOKUP[var]["obs_name"]
        ds[var].attrs["original_long_name"] = obs_LOOKUP[var]["obs_long_name"]

        # rename dimensions
        ds[var] = ds[var].rename({"latitude": "lat", "longitude": "lon"})

        # additional attributes -- hard coded for EOBS
        ds[var].attrs["time_freq"]          = "daily"  # possible values: daily, hourly, monthly, yearly
        ds[var].attrs["spatial_resolution"] = "0.1deg"  
        ds[var].attrs["domain"]             = "europe"

    
    # Soft check for CF compliance 
    cf_status(ds)

    return ds




