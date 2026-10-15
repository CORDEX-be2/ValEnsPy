import xclim
import warnings
import xarray as xr

#Eventually replace this with built in xclim functionality - see https://github.com/Ouranosinc/xclim/issues/2127 and https://github.com/Ouranosinc/xclim/pull/2144
def convert_units_to(data : xr.Dataset | xr.DataTree, var : str, target_unit : str | dict, context: str ="infer"):
    """
    Convert units of a variable in a xr.Dataset or xr.DataTree to a target unit. 
    
    convert_units_to wraps :func:`xclim.units.convert_units_to` to enable usage with xr.Dataset and xr.DataTree objects.

    Parameters
    ----------
    data : xr.Dataset or xr.DataTree
        The input dataset or data tree.
    var : str
        The name of the variable to convert.
    target_unit : str or dict
        The target unit to convert to. A dict ``{"units": ..., "units_metadata": ...}`` is passed through
        to xclim, e.g. with ``"temperature: difference"`` to convert a temperature difference.
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
    
def _ds_convert_units_to(ds: xr.Dataset, var: str, target_unit: str | dict, context: str = "infer"):
    """Convert units of a variable in a Dataset to a target unit."""
    if var in ds:
        ds = ds.copy()
        ds[var] = xclim.units.convert_units_to(ds[var], target_unit, context=context)
    elif ds.data_vars:
        warnings.warn(f"Variable {var} not found in the dataset. Conversion not applied.")
    return ds
    
def _dt_convert_units_to(dt: xr.DataTree, var: str, target_unit: str | dict, context: str = "infer"):
    """Convert units of a variable in a DataTree to a target unit."""
    
    return dt.map_over_datasets(
        _ds_convert_units_to,
        var,
        target_unit,
        context
    )


def convert_units_by_mapping(data : xr.Dataset | xr.DataTree, units_map : dict):
    """
    Convert every variable whose units match a key of `units_map` to that key's target unit.

    Variables are matched on their current ``units`` attribute (exact string), so all variables
    with a given unit are converted without naming them. Variables with other units are untouched.
    For a DataTree, a variable's units are read from the first leaf that has it.

    A variable marked ``units_metadata = "temperature: difference"`` is converted as a temperature
    difference, not an absolute temperature: no offset is applied (K -> degC leaves the values
    unchanged, K -> degF only rescales) and the mark is kept on the result.

    Parameters
    ----------
    data : xr.Dataset or xr.DataTree
        The input dataset or data tree.
    units_map : dict
        ``{source_unit: target}``, where ``target`` is a unit string, a ``(unit, context)`` tuple, or a
        dict ``{"units": ..., "context": ..., "label": ...}`` (``context`` and ``label`` optional), e.g.
        ``{"K": "degC", "kg m-2 s-1": ("mm/d", "hydro"), "K days": {"units": "degC days", "label": "°C days"}}``.
        The ``hydro`` context is needed to relate a mass flux to a depth rate. ``label`` replaces the
        ``units`` attribute after the conversion, for when pint's own formatting (e.g. ``d degC``)
        is not what should be displayed.

    Returns
    -------
    xr.Dataset or xr.DataTree
        A new dataset or data tree, of the same type as the input.

    See Also
    --------
    :func:`convert_units_to`
    """
    if isinstance(data, xr.DataTree):
        datasets = [node.to_dataset() for node in data.subtree if node.has_data]
    elif isinstance(data, xr.Dataset):
        datasets = [data]
    else:
        raise TypeError("Input data must be either an xarray Dataset or a DataTree.")

    attrs_by_var = {}
    for ds in datasets:
        for var, da in ds.data_vars.items():
            attrs_by_var.setdefault(var, da.attrs)

    for var, attrs in attrs_by_var.items():
        target = units_map.get(attrs.get("units"))
        if target is None:
            continue
        label = None
        if isinstance(target, str):
            target_unit, context = target, "infer"
        elif isinstance(target, dict):
            target_unit, context, label = target["units"], target.get("context", "infer"), target.get("label")
        else:
            target_unit, context = target
        if attrs.get("units_metadata") == "temperature: difference":
            target_unit = {"units": target_unit, "units_metadata": "temperature: difference"}
        data = convert_units_to(data, var=var, target_unit=target_unit, context=context)
        if label is not None:
            data = _set_units_label(data, var, label)
    return data


def _set_units_label(data, var, label):
    """Set the ``units`` attribute of `var` to `label` in a Dataset or every leaf of a DataTree."""
    def _ds(ds):
        return ds.assign({var: ds[var].assign_attrs(units=label)}) if var in ds else ds

    return data.map_over_datasets(_ds) if isinstance(data, xr.DataTree) else _ds(data)
