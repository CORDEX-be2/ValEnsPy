# Create functionality which checks if the an xarray dataset is a ValEnsPy CF compliant dataset
import xarray as xr
from typing import Union, List
from valenspy._utilities import load_yml, load_xarray_from_data_sources
from pathlib import Path

CORDEX_VARIABLES = load_yml("CORDEX_variables")

# Expected metadata attributes
MAIN_METADATA = ["Conventions", "history"]
VARIABLE_METADATA = ["units", "standard_name", "long_name"]


def cf_status(netCDF: Union[str, Path, xr.Dataset]) -> None:
    """
    Provides an overview of the degree of CF compliance of the netCDF file.

    Parameters
    ----------
    netCDF : Union[str, Path, xr.Dataset]
        The netCDF file to check or the xarray dataset to check
    """

    ds = load_xarray_from_data_sources(netCDF)

    print(
        "The file is {status}ValEnsPy CF compliant.".format(
            status="NOT " if not is_cf_compliant(ds, verbose=True) else ""
        )
    )

    non_cf_compliant_vars = []
    cf_compliant_vars = []
    unknown_vars = []
    for var in ds.data_vars:
        if _check_variable_by_name(ds[var]):
            cf_compliant_vars.append(var)
        else:
            if var in CORDEX_VARIABLES:
                non_cf_compliant_vars.append(var)
            else:
                unknown_vars.append(var)

    print(
        f"{len(cf_compliant_vars)/len(ds.data_vars)*100:.2f}% of the variables are ValEnsPy CF compliant"
    )
    if cf_compliant_vars:
        print(f"ValEnsPy CF compliant: {cf_compliant_vars}")
    if non_cf_compliant_vars:
        print(f"NOT ValEnsPy CF compliant: {non_cf_compliant_vars}")
    if unknown_vars:
        print(f"Unknown to ValEnsPy: {unknown_vars}")

    for var in non_cf_compliant_vars:
        print(
            f"The following attributes are missing or incorrect for the variable {var}:"
        )
        print("{:<15} {:<25} {:<25}".format("Attribute", "Actual", "Expected"))
        print("-" * 65)
        for attr in VARIABLE_METADATA:
            actual = ds[var].attrs.get(attr, "Not present")
            expected = CORDEX_VARIABLES.get(var).get(attr, "Not present")
            if actual != expected:
                print("{:<15} {:<25} {:<25}".format(attr, actual, expected))


def is_cf_compliant(netCDF: Union[str, Path, xr.Dataset], verbose=False) -> bool:
    """
    Check if a file is a ValEnsPy CF compliant netCDF file. The following checks are performed:
    - Check if the file is a netCDF file (or an xarray dataset)
    - Check if the main metadata attributes exist (title, history)
    - Check if the variable metadata attributes exist (units, standard_name, long_name)
    - For each variable check if it is present in the predefined variables, if so check if the attributes match.

    Parameters
    ----------
    netCDF : Union[str, Path, xr.Dataset]
        The netCDF file to check or the xarray dataset to check
    verbose : bool, optional
        If True, print the results of the checks, by default False

    Returns
    -------
    bool
        True if the file is ValEnsPy CF compliant, False otherwise

    Examples
    --------
    .. code-block:: python

        >>> import valenspy as vp
        >>>
        >>> vp.is_cf_compliant(vp.demo_data_CF)
        True

    """
    ds = load_xarray_from_data_sources(netCDF)

    var_meta_data_ok = all(
        [
            _check_variable_metadata(ds[var])
            for var in ds.data_vars
            if "_bnds" not in var
        ]
    )
    main_meta_data_ok = _check_main_metadata(ds)
    cordex_vars_data_ok = all(
        [
            _check_variable_by_name(ds[var])
            for var in ds.data_vars
            if var in CORDEX_VARIABLES
        ]
    )
    time_dimension_ok = _check_time_dimension(ds)

    if verbose:
        if not var_meta_data_ok:
            print("Variable metadata is missing or incorrect")
        if not main_meta_data_ok:
            print("Main metadata is missing or incorrect")
        if not cordex_vars_data_ok:
            print("Variable attributes are missing or incorrect")
        if not time_dimension_ok:
            print("Time dimension is missing or has an incorrect type")

    return (
        var_meta_data_ok
        and main_meta_data_ok
        and cordex_vars_data_ok
        and time_dimension_ok
    )


def _check_time_dimension(ds: xr.Dataset):
    """
    Check if the time dimension is present and is of the type datetime64

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to check

    Returns
    -------
    bool
        True if the time dimension is present and is of the type datetime64, False otherwise
    """
    return "time" in ds.dims and ds.time.dtype == "datetime64[ns]"


def _check_variable_by_name(da: xr.DataArray):
    """
    Check if the variable is a CORDEX variable.
    If it is check if the attributes in the datarray match the corresponding attributes in the CORDEX_VARIABLES dictionary

    Parameters
    ----------
    da : xr.DataArray
        The data array to check

    Returns
    -------
    bool
        True if the variable is a CORDEX variable and the attributes match or if the variable is not a CORDEX variable, False otherwise
    """
    if da.name in CORDEX_VARIABLES:
        return all(
            [
                CORDEX_VARIABLES[da.name][attr] == da.attrs[attr]
                for attr in da.attrs
                if attr in CORDEX_VARIABLES[da.name]
            ]
        )
    else:
        return False


def _check_main_metadata(ds: xr.Dataset):
    """
    Check if the required main metadata attributes exist

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to check

    Returns
    -------
    bool
        True if the required main metadata attributes exist, False otherwise
    """

    return all([attr in ds.attrs for attr in MAIN_METADATA])


def _check_variable_metadata(da: xr.DataArray):
    """
    Check if the CF required variable attributes exist

    Parameters
    ----------
    da : xr.DataArray
        The data array to check

    Returns
    -------
    bool
        True if the CF required variable attributes exist, False otherwise
    """

    return all([attr in da.attrs for attr in VARIABLE_METADATA])
