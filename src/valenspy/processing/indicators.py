from xarray import DataTree
import xarray as xr
from xclim.core.indicator import Indicator

#Eventually replace this with built in xclim functionality - see https://github.com/Ouranosinc/xclim/issues/2127 and https://github.com/Ouranosinc/xclim/pull/2144
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
    def xclim_indicator_ds(ds: xr.Dataset, indicator: Indicator,  vars: str | list, **kwargs) -> xr.Dataset:
        """
        Wrapper function to apply the indicator to a single dataset.
        """
        if isinstance(vars, str):
            if vars in ds:
                return indicator(ds[vars], **kwargs).to_dataset()
        elif isinstance(vars, list): #Order is important!
            if all(var in ds for var in vars):
                data_arrays = [ds[var] for var in vars]
                return indicator(*data_arrays, **kwargs).to_dataset()
        
    return dt.map_over_datasets(
        xclim_indicator_ds,
        indicator,
        vars,
        kwargs=kwargs
    )
