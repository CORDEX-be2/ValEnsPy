import regionmask
import xarray as xr

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


