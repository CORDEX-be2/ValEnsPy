import regionmask
import xarray as xr
import geopandas as gpd
from pathlib import Path

def add_regionmask_region(ds, region_mask_region):
    """
    Add a region from the regionmask package to a dataset. Regions will be added as a dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to add region to.
    region_mask_region : regionmask.Regions
        Region from the regionmask package. e.g. regionmask.defined_regions.prudence

    Returns
    -------
    xarray.Dataset
        Dataset with region as a new dimension.
    """
    
    mask = region_mask_region.mask_3D(ds.lon, ds.lat)
    return ds.where(mask)

def create_regionmask_region(shapefile_path, abbrevs=None, name=None):
    """
    Create a region from a shapefile using the regionmask package.

    Parameters
    ----------
    shapefile_path : Path
        Path to the shapefile.
    abbrevs : str, optional
        The column name in the shapefile that contains the region abbreviations.
    name : str, optional
        The column name in the shapefile that contains the region names.

    Returns
    -------
    regionmask.Regions
        Region from the regionmask package.
    """
    
    gdf_shp = gpd.read_file(shapefile_path)
    return regionmask.from_geopandas(gdf_shp, abbrevs=abbrevs, name=name)

def get_shapefile_mask(ds: xr.Dataset, shapefile_path: Path):

    """
    Generates a mask from a shapefile to apply to an xarray Dataset.

    This function reads a shapefile using Geopandas, converts it to the WGS84 coordinate reference system (CRS),
    and creates a mask that can be applied to the input xarray Dataset. The mask identifies the grid cells that
    fall within the shapefile's region.

    Parameters
    ----------
    ds : xr.Dataset
        Input xarray Dataset containing longitude and latitude coordinates.
    shapefile_path : Path
        Path to the shapefile to be used for masking.

    Returns
    -------
    mask_shp : xr.Dataset
        A boolean mask array where grid cells within the shapefile region are marked as True,
        and those outside are marked as False.

    Notes
    -----
    - The shapefile is converted to the WGS84 CRS (EPSG:4326) before creating the mask.
    - The function uses the regionmask library to generate the mask.

    Examples
    --------
    >>> import xarray as xr
    >>> from pathlib import Path
    >>> ds = xr.open_dataset('path_to_your_dataset.nc')
    >>> shapefile = Path('path_to_your_shapefile.shp')
    >>> mask = get_shapefile_mask(ds, shapefile)
    """
    # read shape file into geopandas geodataframe
    gdf_shp = gpd.read_file(shapefile_path)

    # convert geodataframe to WGS84 and mask - only needed to do once
    gdf_shp = gdf_shp.to_crs(epsg=4326)

    # do masking
    mask_shp = (
        regionmask.mask_geopandas(gdf_shp, ds.lon.values, ds.lat.values) + 1
    ) > 0

    return mask_shp
