"""A collection of functions unique to an input dataset required for conversion to valenspy compliant xarray datasets."""

from valenspy._utilities import load_yml, _fix_lat_lon
import xarray as xr
import xarray.coding.times as _xr_times

CORDEX_VARIABLES = load_yml("CORDEX_variables")


def _fix_noresm2mm_calendar(ds: xr.Dataset) -> xr.Dataset:
    """MAR's NorESM2-MM-driven raw files (source: .../MARv3.14-NorESM2-MM-...) are
    internally 365-day/noleap data (confirmed against the raw files directly: a leap
    year has exactly 365 timesteps, vs. 366 for MAR's other driving models,
    EC-Earth3-Veg and MPI-ESM1-2-HR, in the same year) but declare calendar=
    "standard" in their metadata. Decoded normally, that silently drops Feb 29 as a
    real gap rather than a calendar-defined absence, which breaks any downstream
    frequency check (e.g. xclim's indicator datachecks) that spans a leap year.

    Re-derives the correct dates by re-encoding the (wrongly-decoded) time values
    back to the file's own true raw numeric offsets - using the exact units/calendar
    xarray recorded in .encoding while decoding, not an auto-chosen reference; passing
    units=None to encode_cf_datetime re-anchors to the array's own first value and
    silently turns this into a no-op, since any single already-decoded file's own
    timestamps are always internally regular by themselves - the drift only exists
    relative to the *true* calendar, so recovering it requires the file's real
    original reference date, which only survives in .encoding, not in the values
    array. Then decodes those same raw numbers again under the correct "noleap"
    calendar - verified numerically identical to decoding straight from the raw file
    with the calendar attribute corrected before decoding, for every year checked.

    Only valid applied per file, before any cross-file concatenation: MAR's
    NorESM2-MM files are delivered in ~5-year batches, each with its own "hours
    since <epoch>" reference that resets rather than continuing from one fixed
    origin, so the same correction applied to an already-concatenated, multi-file
    series can't tell which portion came from which batch and silently produces
    wrong dates - this must run here, in the per-file Input Convertor, not later in
    any per-ensemble post-processing.

    Deliberately kept as cftime.DatetimeNoLeap (not converted to datetime64[ns]):
    a noleap date is always a valid Gregorian date too, so datetime64 could
    represent the values, but not the calendar - a downstream frequency check
    (e.g. xclim's) needs the calendar tag itself to recognize "no Feb 29" as
    expected rather than a gap; plain datetime64 values with Feb 29 skipped look
    just as irregular to such a check as the original bug did. See
    _convert_all_units_to_CF in _utilities/unit_converter.py, which used to force
    everything through pd.to_datetime() (breaking on cftime) - fixed alongside this
    to skip that conversion when the time coordinate is already cftime-indexed.
    """
    if "time" not in ds.coords:
        return ds
    source = ds.encoding.get("source", "")
    if "NorESM2-MM" not in source:
        return ds
    # Must reuse the file's own true units/calendar from .encoding (preserved by
    # xarray's decoder even after decoding), not an auto-chosen reference - passing
    # units=None re-anchors to the array's own first value and silently turns this
    # into a no-op, since a single already-decoded file's own timestamps are always
    # internally regular (the drift only shows up relative to the *true* calendar,
    # never within the wrongly-decoded values by themselves).
    time_var = ds["time"]
    raw_numeric, units, _ = _xr_times.encode_cf_datetime(
        time_var.values, units=time_var.encoding["units"], calendar=time_var.encoding["calendar"]
    )
    corrected_time = _xr_times.decode_cf_datetime(raw_numeric, units=units, calendar="noleap")
    return ds.assign_coords(time=corrected_time)


def EOBS_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw EOBS data to a ValensPy compliant xarray Dataset.

    Rename latitude and longitude coordinates to lat and lon, respectively.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of EOBS observations to convert

    Returns
    -------
    Dataset
        The CF compliant EOBS observations for the specified variable.
    """
    ds = _fix_lat_lon(ds)

    return ds


def ERA5_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw ERA5 data to a ValensPy compliant xarray Dataset.

    Rename latitude and longitude coordinates to lat and lon, respectively. Rename valid_time to time for certain variables.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ERA5 observations to convert

    Returns
    -------
    Dataset
        The CF compliant ERA5 observations for the specified variable.
    """

    # bugfix ERA5 (found in clh): replace valid_time by time
    if "time" not in ds:
        ds = ds.rename({"valid_time": "time"})

    ds = _fix_lat_lon(ds)

    return ds


def CCLM_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw CCLM data to a ValensPy compliant xarray Dataset.

    Flatten the pressure dimension by renaming variables with the pressure level in the name and removing the pressure dimension.
    Drop the last time step of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CCLM simulation to convert

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """

    #For each variable in the dataset which has a pressure dimension, create a new variable with the pressure level in the name and remove the pressure dimension
    if "pressure" in ds.dims:
        for var in ds.data_vars:
            if "pressure" in ds[var].dims:
                for pressure in ds[var].pressure.values:
                    new_var = var + str(int(pressure / 100)) + "p"
                    ds[new_var] = ds[var].sel(pressure=pressure)
                ds = ds.drop_vars(var)
        ds = ds.drop_dims("pressure")
        
    return ds


def ALARO_K_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw ALARO_K data to a ValensPy compliant xarray Dataset.

    Does nothing, WIP

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of ALARO_K simulation to convert

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """
    # # Special conversion for precipitation
    # if "rain_convective" in ds.data_vars and "rain_stratiform" in ds.data_vars:
    #     ds["pr"] = ds["rain_convective"] + ds["rain_stratiform"]
    #     ds["pr"].attrs["units"] = "mm"
    #     ds["pr"] = convert_mm_to_kg_m2s(ds["rain_convective"] + ds["rain_stratiform"])
    #     ds["pr"].attrs["standard_name"] = "precipitation_flux"
    #     ds["pr"].attrs["long_name"] = "Precipitation"
    #     ds["pr"].attrs["dataset"] = model_name
    #     ds["pr"].attrs["original_name"] = "rain_convective + rain_stratiform"
    #     for key, value in metadata_info.items():
    #         ds["pr"].attrs[key] = value

    #     # Assuming monthly decumilation! This is not always the case!
    #     def decumilate(ds):
    #         ds_decum = ds.diff("time")
    #         # Add the first value of the month of original dataset to the decumilated dataset
    #         ds_decum = xr.concat([ds.isel(time=0), ds_decum], dim="time")
    #         return ds_decum

    #     ds.coords["year_month"] = ds["time.year"] * 100 + ds["time.month"]
    #     ds["pr"] = ds["pr"].groupby("year_month").apply(decumilate)

    return ds

def RADCLIM_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw RADCLIM data to a ValensPy compliant xarray Dataset.

    Rename nlon and nlat to lon and lat, respectively. Set the coordinates lat_bounds and lon_bounds as coordinates.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of CCLM simulation to convert

    Returns
    -------
    Dataset
        The CF compliant CCLM model data for the specified variable.
    """

    ds = ds.set_coords(("lat_bounds", "lon_bounds"))

    if "nlon" in ds.dims:
        ds = ds.rename({"nlon": "lon"})
    if "nlat" in ds.dims:
        ds = ds.rename({"nlat": "lat"})

    return ds

def MAR_to_CF(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert a xarray with raw MAR data to a ValensPy compliant xarray Dataset.

    Rename TIME to time and remove the ZTQLEV and ZUVLEV dimensions by selecting the first value of each dimension.
    Derive scalar near-surface wind speed (sfcWind) from the U2Z/V2Z wind components, since MAR - unlike every
    other source in INPUT_CONVERTORS - does not report a scalar wind speed variable directly.
    For files driven by NorESM2-MM specifically, also corrects a mislabeled calendar
    (declared "standard", actually 365-day/noleap) - see _fix_noresm2mm_calendar.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset of MAR simulation to convert

    Returns
    -------
    Dataset
        xarray dataset ready for unit conversion

    """
    ds = ds.rename({'TIME':'time'})
    ds = _fix_noresm2mm_calendar(ds)
    ds = ds.isel(ZTQLEV=0,ZUVLEV=0)

    # MAR only provides the wind vector components (U2Z, V2Z), not a scalar
    # wind speed - derive it here so it can go through the same raw_name ->
    # CORDEX rename/unit-conversion path as every other variable (see the
    # "sfcWind" entry in MAR_lookup.yml). Computed before the raw_name ->
    # CORDEX rename, at U2Z/V2Z's native "m/s", matching their raw_units in
    # the lookup table.
    if "U2Z" in ds and "V2Z" in ds:
        ds["sfcWind_derived"] = (ds["U2Z"] ** 2 + ds["V2Z"] ** 2) ** 0.5
        ds["sfcWind_derived"].attrs["units"] = "m/s"

    return ds