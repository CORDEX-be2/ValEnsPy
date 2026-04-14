import regionmask
import xarray as xr

def mask_dt(dt : xr.DataTree, mask : xr.Dataset | xr.DataArray) -> xr.DataTree:
    """
    Mask a datatree with a given mask. The mask will be applied to all datasets in the datatree.

    Parameters
    ----------
    dt : xarray.DataTree
        DataTree to be masked.
    mask : xarray.Dataset or xarray.DataArray
        Mask to be applied to the datatree. The mask should have the same dimensions as the datasets in the datatree.

    Returns
    -------
    xarray.DataTree
        Masked datatree.
    """
    return dt.map_over_datasets(lambda ds: ds.where(mask) if ds else ds) #Trick to deal with leaves in the datatree that are not datasets (e.g. groups)

def add_prudence_regions(ds : xr.Dataset) -> xr.Dataset:
    """
    Add PRUDENCE regions to a dataset. Regions will be added as a dimension (3D mask).
    The PRUDENCE regions are defined in the regionmask package.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to add PRUDENCE regions to.

    Returns
    -------
    xarray.Dataset
        Dataset with PRUDENCE regions as a new dimension.
    """
    
    prudence = regionmask.defined_regions.prudence
    mask = prudence.mask_3D(ds.lon, ds.lat)
    return ds.where(mask)


