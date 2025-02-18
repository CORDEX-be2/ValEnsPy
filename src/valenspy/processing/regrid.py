import cdo
import xarray as xr
import cf_xarray
import xesmf as xe

def remap_xesmf(ds, ds_out, method="bilinear", regridder_kwargs={}, regridding_kwargs={}):
    """Remap the input dataset to the target grid using xESMF.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset to remap.
    ds_out : xarray.Dataset
        The target grid dataset.
    method : str, optional
        The remap method to use, by default "bilinear".
    regridder_kwargs : dict
        Keyword arguments for the creation of the regridder. See xe.Regridder for more information.
    regridding_kwargs : dict
        Keyword arguments for the actual regridding. See xe.Regridder.regrid_dataset() for more information.

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

def remap_cdo(ds, target_grid, remap_method="bil", tempdir=None, output_path=None):
    """Remap the input dataset to the target grid using CDO.

    Parameters
    ----------
    target_grid : str
        The target grid file. This is either in the form "r360x180" for a regular global grid longitude x latitude,
        a path to a netCDF file with the target grid or a path to a text file with the target grid.
    ds : xarray.Dataset
        The input dataset to remap. Note in the background this is saved to disk as a netCDF file, processed by CDO and read back in as an xarray.Dataset.
        This can be memory intensive for large datasets.
    remap_method : str
        The remap method to use.
    tempdir : str, optional
        The temporary directory to save the temporary files, by default None (uses /tmp).
    output_path : bool, optional
        If True, the output file is saved to disk, by default False.

    Returns
    -------
    xarray.Dataset
        The remapped dataset in xarray format.
    """
    if tempdir:
        cdo_obj = cdo.Cdo(tempdir=tempdir)
    else:
        cdo_obj = cdo.Cdo()

    if remap_method == "bil":
        remap = cdo_obj.remapbil(target_grid, input=ds, returnXDataset=True)
    elif remap_method == "con":
        # for conservative remapping, lat_bounds and lon_bounds are required. Check whether these are present and if not add these.
        if not ("lat_bounds" in ds.variables and "lon_bounds" in ds.variables):
            ds = ds.cf.add_bounds(("lat", "lon"))
        remap = cdo_obj.remapcon(target_grid, input=ds, returnXDataset=True)
    elif remap_method == "dis":
        remap = cdo_obj.remapdis(target_grid, input=ds, returnXDataset=True)
    elif remap_method == "nn":
        remap = cdo_obj.remapnn(target_grid, input=ds, returnXDataset=True)
    else:
        raise ValueError(
            f"Remap method {remap_method} not supported: choose from 'bil', 'con', 'dis' or 'nn'."
        )
    if output_path:
        remap.to_netcdf(output_path)
    # Make sure the remap dataset is an dask array
    remap = remap.chunk(chunks="auto")

    # Rename the longitude and latitude to lon and lat if they are not already named like this.
    if "longitude" in remap:
        remap = remap.rename({"longitude": "lon"})
    if "latitude" in remap:
        remap = remap.rename({"latitude": "lat"})

    # make sure lat and lon are sorted ascending
    remap = remap.sortby("lat").sortby("lon")

    return remap
