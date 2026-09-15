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


def _restrict_to_shared_vars(ds: xr.Dataset, other_vars) -> xr.Dataset:
    """`ds` restricted to whichever of its own variables are also in `other_vars`."""
    return ds[[v for v in ds.data_vars if v in other_vars]]


def _mask_ds_to_reference_coverage(ds: xr.Dataset, reference: xr.Dataset, dim) -> xr.Dataset:
    """One dataset's worth of `mask_to_reference_coverage` - used directly by
    `Model2Self.apply` (a single Dataset, no `mask_dt` to delegate to there)."""
    coverage = reference.notnull().any(dim=dim)
    return _restrict_to_shared_vars(ds, reference.data_vars).where(coverage)


def mask_to_reference_coverage(dt: xr.DataTree, reference: xr.Dataset, dim: str = "time") -> xr.DataTree:
    """Mask every leaf of `dt` to `reference`'s own coverage: NaN wherever `reference`
    has no valid value anywhere along `dim` (default "time", i.e. mask by spatial
    extent - pass e.g. `["lat", "lon"]` for temporal coverage instead).

    Each leaf is first restricted to whichever of its own variables `reference` also
    has - a variable `reference` doesn't cover at all is dropped rather than silently
    masked or left untouched (plain `.where(reference)` would otherwise drop it anyway,
    with no error) - so this works whether `dt`/`reference` carry one variable or
    several (e.g. AnnualCycle, applied to every variable at once). Restrict `reference`
    to one variable first (e.g. `reference[[var]]`) to mask only that variable.

    Returns
    -------
    xr.DataTree
        `dt`, masked to `reference`'s own coverage.
    """
    coverage = reference.notnull().any(dim=dim)
    restricted = dt.map_over_datasets(
        lambda ds: _restrict_to_shared_vars(ds, reference.data_vars) if ds else ds
    )
    return mask_dt(restricted, coverage)

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


