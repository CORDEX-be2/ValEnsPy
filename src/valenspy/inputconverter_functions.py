from pathlib import Path
import xarray as xr


def _non_convertor(file: Path) -> Path:
    """A dummy function that does not convert the input file as it already is in the correct format."""
    return file



def EOBS_to_CF(file: Path) -> Path:
    """Convert the EOBS netCDF file to a netCDF file in CF convention, with all updated coordinate names, attributes and 
        unit conversions."""

    ds = xr.open_mfdataset(file, combine='by_coords', chunks='auto')

    # make EOBS CF compliant (and CMORized)
    ds = ds.rename_vars({EOBS_variable: variable})
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})

    # Convert from Celsius to Kelvin -- to do: put this in a seperate convertor to be used by other functions? 
    # or create a look-up table for observations including units, and based thereon necessary conversions to be done? 
    if variable == 'tas': 
        ds[variable] = ds[variable] + 273.15

    # add necessary metadata
    ds.attrs["long_name"]          = df_cf_lookup.loc[df_cf_lookup['variable_name'] == variable,'long_name'].values[0] # from the CORDEX look-up table
    ds.attrs["original_name"]      = EOBS_variable
    ds.attrs["original_long_name"] = df_lookup_EOBS.loc[df_lookup_EOBS['variable_name'] == variable,'EOBS_long_name'].values[0]
    ds.attrs["units"]              = df_cf_lookup.loc[df_cf_lookup['variable_name'] == variable,'units'].values[0] # from the CORDEX look-up table 
    ds.attrs["time_freq"]          = "daily"  # possible values: daily, hourly, monthly, yearly
    ds.attrs["spatial_resolution"] = "0.1deg"  
    ds.attrs["domain"]             = "europe"

    return ds
