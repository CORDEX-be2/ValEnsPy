from datatree import map_over_subtree, DataTree
from xclim.core.indicator import Indicator

#Note - map_over_subtree is from the original xarray-datatree package, when moved to xarray it was not kept. When shifting to datatree from xarray this should be dealt with.
@map_over_subtree
def xclim_indicator(dt : DataTree, indicator : Indicator, vars : str | list, **kwargs) -> DataTree:
    """
    Calculate an xclim indicator on a data tree.

    This function is a wrapper around the xclim indicator functions. It takes a data tree and applies the indicator to each dataset in the tree.
    The indicator function should be a :py:class:`xclim.indicator` and the variables should match the expected order of that indicator function.

    Parameters
    ----------
    dt : DataTree
        A data tree containing the datasets to calculate the indicator on.
    indicator : :py:class:`xclim.indicator`
        An xclim indicator function.
    vars : str or list
        The variable(s) to calculate the indicator with. The order of the variables is important and should match the order expected by the xclim indicator.
    **kwargs
        Additional keyword arguments to pass to the indicator function.

    Returns
    -------
    xr.Dataset
        A new dataset with the indicator calculated
    """
    #Note dt is actually a dataset here but due to the wrapper the function works on the datatree level by broadcasting this function to all datasets in the tree.
    if isinstance(vars, str):
        return indicator(dt[vars], **kwargs).to_dataset()
    elif isinstance(vars, list): #Order is important!
        data_arrays = [dt[var] for var in vars]
        return indicator(*data_arrays, **kwargs).to_dataset()