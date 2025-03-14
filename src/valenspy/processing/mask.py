def add_prudence_regions(ds):
    """
    Add PRUDENCE regions to a dataset. Regions will be added as a dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to add PRUDENCE regions to.

    Returns
    -------
    xarray.Dataset
        Dataset with PRUDENCE regions as a new dimension.
    """

    import regionmask
    prudence = regionmask.defined_regions.prudence
    mask = prudence.mask_3D(ds.lon, ds.lat)
    return ds.where(mask)


