from pathlib import Path
import xarray as xr


def _non_convertor(file: Path) -> Path:
    """A dummy function that does not convert the input file as it already is in the correct format."""
    return file


# work in progress. 
def EOBS_to_CF(directory: Path, variable: string) -> Path:
    """Convert the EOBS netCDF file to a netCDF file in CF convention for specified CF variable."""

    # load look-up table of EOBS variables to CORDEX variables
    df_lookup_EOBS = pd.read_csv('../ancilliary_data/lookuptable_EOBS.csv', sep=";")

    # get E-OBS variable corresponding to the requested variable using its look-up table
    EOBS_variable = df_lookup_EOBS.loc[df_lookup_EOBS['variable_name'] == variable,'EOBS_name'].values[0]

    # open the EOBS file for the corresponding variable
    files = list(directory.glob("*"+EOBS_variable+"*mean*.nc")) #Select all the netCDF files in the directory

    ds = xr.open_mfdataset(files, combine='by_coords', chunks='auto')

    # make EOBS CF compliant (and CMORized)
    ds = ds.rename_vars({EOBS_variable: variable})
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})

    # Convert from Celsius to Kelvin -- to do: put this in a seperate convertor to be used by other functions? 
    if variable == 'tas': 
        ds[variable] = ds[variable] + 273.15

    # add necessary metadata -- this makes the person including the obs data think whether additional processing is needed 

    ds.attrs["long_name"] = df_cf_lookup.loc[df_cf_lookup['variable_name'] == variable,'long_name'].values[0] # from the CORDEX look-up table
    ds.attrs["original_name"] = EOBS_variable
    ds.attrs["original_long_name"] = df_lookup_EOBS.loc[df_lookup_EOBS['variable_name'] == variable,'EOBS_long_name'].values[0]
    ds.attrs["units"] = df_cf_lookup.loc[df_cf_lookup['variable_name'] == variable,'units'].values[0] # from the CORDEX look-up table 
    ds.attrs["time_freq"] = "daily"  # possible values: daily, hourly, monthly, yearly
    ds.attrs["spatial_resolution"] = "0.1deg"  
    ds.attrs["domain"] = "europe"

    return ds
