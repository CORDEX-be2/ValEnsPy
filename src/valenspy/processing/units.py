import xclim
import warnings
import xarray as xr

#Eventually replace this with built in xclim functionality - see https://github.com/Ouranosinc/xclim/issues/2127 and https://github.com/Ouranosinc/xclim/pull/2144
def convert_units_to(data : xr.Dataset | xr.DataTree, var : str, target_unit : str, context: str ="infer"):
    """
    Convert units of a variable in a xr.Dataset or xr.DataTree to a target unit. 
    
    convert_units_to wraps :func:`xclim.units.convert_units_to` to enable usage with xr.Dataset and xr.DataTree objects.

    Parameters
    ----------
    data : xr.Dataset or xr.DataTree
        The input dataset or data tree.
    var : str
        The name of the variable to convert.
    target_unit : str
        The target unit to convert to.
    context : str, optional
        The context in which the conversion is made. Default is 'infer'.
    
    Returns
    -------
    xr.Dataset or xr.DataTree
        A new dataset or data tree with the variable converted to the target unit depending on the input type.

    See Also
    --------
    :func:`xclim.units.convert_units_to`
    """

    if isinstance(data, xr.DataTree):
        return _dt_convert_units_to(data, var, target_unit, context)
    elif isinstance(data, xr.Dataset):
        return _ds_convert_units_to(data, var, target_unit, context)
    else:
        raise TypeError("Input data must be either an xarray Dataset or a DataTree.")
    
def _ds_convert_units_to(ds: xr.Dataset, var: str, target_unit: str, context: str = "infer"):
    """Convert units of a variable in a Dataset to a target unit."""
    if var in ds:
        ds = ds.copy()
        ds[var] = xclim.units.convert_units_to(ds[var], target_unit, context=context)
    elif ds.data_vars:
        warnings.warn(f"Variable {var} not found in the dataset. Conversion not applied.")
    return ds
    
def _dt_convert_units_to(dt: xr.DataTree, var: str, target_unit: str, context: str = "infer"):
    """Convert units of a variable in a DataTree to a target unit."""
    
    return dt.map_over_datasets(
        _ds_convert_units_to,
        var,
        target_unit,
        context
    )
