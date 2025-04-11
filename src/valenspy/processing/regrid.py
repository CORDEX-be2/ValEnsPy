import xarray as xr
import cf_xarray
import xesmf as xe

def remap_xesmf(ds : xr.Dataset, ds_out : xr.Dataset, method : str ="bilinear", regridder_kwargs : dict ={}, regridding_kwargs: dict ={}):
    """Remap the input dataset to the target grid using xESMF.

    If lat_bounds and lon_bounds are not present in the input dataset, they will be added.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset to remap.
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
    xarray.Dataset
        The remapped dataset in xarray format.
    """
    if method=="conservative":
        if not ("lat_bounds" in ds.variables and "lon_bounds" in ds.variables):
            ds = ds.cf.add_bounds(("lat", "lon"))
    regridder = xe.Regridder(ds, ds_out, method, **regridder_kwargs)
    ds_reg = regridder(ds, **regridding_kwargs)
    return ds_reg
