import xclim
import warnings
import xarray as xr
from datatree import DataTree, map_over_subtree

def convert_units_to(data : xr.Dataset | DataTree, var : str, target_unit : str, context: str ="infer"):
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
    if isinstance(data, DataTree):
        return _convert_units_dt_to(data, var, target_unit, context)
    elif isinstance(data, xr.Dataset) and var in data:
        ds = data.copy()
        ds[var] = xclim.units.convert_units_to(ds[var], target_unit, context=context)
        return ds
    else:
        warnings.warn(f"Variable {var} not found in the dataset or data tree.")
        return data

@map_over_subtree
def _convert_units_dt_to(ds, var, target_unit, context="infer"):
    if var in ds:
        ds = ds.copy()
        ds[var] = xclim.units.convert_units_to(ds[var], target_unit, context=context)
    return ds
