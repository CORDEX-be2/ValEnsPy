import xarray as xr
import cf_xarray
try:
    # Import check for xesmf as it has non-python dependencies. This ensures valenspy can function without xesmf - for example when installed via pip.
    import xesmf as xe
    XESM_AVAILABLE = True
except ImportError:
    XESM_AVAILABLE = False

def remap_xesmf(data : xr.Dataset | xr.DataTree, ds_out : xr.Dataset, method : str ="bilinear", regridder_kwargs : dict ={}, regridding_kwargs: dict ={}):
    """Remap the input dataset to the target grid using xESMF.

    If lat_bounds and lon_bounds are not present in the input dataset, they will be added.

    Parameters
    ----------
    data : xarray.Dataset or xarray.DataTree
        The input data to remap.
    ds_out : xarray.Dataset
        The target grid dataset.
    method : str, optional
        The remap method to use, by default "bilinear".
    regridder_kwargs : dict
        Keyword arguments for the creation of the regridder. See :func:`xesmf.Regridder`.
    regridding_kwargs : dict
        Keyword arguments for the actual regridding. See :func:`xesmf.Regridder.regrid_dataset()`.

    Returns
    -------
    xarray.Dataset or xarray.DataTree
        The remapped data.
    """
    if not XESM_AVAILABLE:
        raise ImportError(f"The xesmf dependency ESMF and esmpy (EMSF's python interface) are not installed. Please install it with 'conda install -c conda-forge esmpy' or similar to use this function.")
    if isinstance(data, xr.DataTree):
        return _remap_xesmf_dt(data, ds_out, method, regridder_kwargs, regridding_kwargs)
    elif isinstance(data, xr.Dataset):
        return _remap_xesmf_ds(data, ds_out, method, regridder_kwargs, regridding_kwargs)
    else:
        raise TypeError("Input data must be either an xarray Dataset or a DataTree.")

def _remap_xesmf_dt(dt : xr.DataTree, ds_out : xr.Dataset, method : str ="bilinear", regridder_kwargs : dict ={}, regridding_kwargs: dict ={}):
    """
    Remap the input DataTree to the target grid using xESMF. Applies _remap_xesmf_ds to each dataset in the DataTree.
    """

    return dt.map_over_datasets(
        _remap_xesmf_ds,
        ds_out,
        method,
        regridder_kwargs,
        regridding_kwargs
    )

def _remap_xesmf_ds(ds : xr.Dataset, ds_out : xr.Dataset, method : str ="bilinear", regridder_kwargs : dict ={}, regridding_kwargs: dict ={}):
    """Remap the input dataset to the target grid using xESMF.

    If lat_bounds and lon_bounds are not present in the input dataset, they will be added.
    """
    if not ds.data_vars:
        return ds
    if method=="conservative":
        if not ("lat_bounds" in ds.variables and "lon_bounds" in ds.variables):
            ds = ds.cf.add_bounds(("lat", "lon"))
    regridder = xe.Regridder(ds, ds_out, method, **regridder_kwargs)
    ds_reg = regridder(ds, **regridding_kwargs)
    return ds_reg
