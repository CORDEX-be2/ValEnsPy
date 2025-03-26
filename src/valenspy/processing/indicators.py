from datatree import map_over_subtree

@map_over_subtree
def xclim_indicator(ds, indicator, vars, **kwargs):
    """
    Calculate an xclim indicator on a data tree.

    Parameters
    ----------
    ds : xr.Dataset
        The input dataset.
    indicator : xclim.indicator
        An xclim indicator function.
    vars : str or list
        The variable(s) to calculate the indicator with. The order of the variables is important.
    **kwargs
        Additional keyword arguments to pass to the indicator function.

    Returns
    -------
    xr.Dataset
        A new dataset with the indicator calculated
    """

    if isinstance(vars, str):
        return indicator(ds[vars], **kwargs).to_dataset()
    elif isinstance(vars, list): #Order is important!
        data_arrays = [ds[var] for var in vars]
        return indicator(*data_arrays, **kwargs).to_dataset()