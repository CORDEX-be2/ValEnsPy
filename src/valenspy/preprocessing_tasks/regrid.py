import cdo
import xarray as xr

def remap_cdo(target_grid, ds, remap_method="bil", output_path=None):
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
    output_file : bool, optional
        If True, the output file is saved to disk, by default False.
    
    Returns
    -------
    xarray.Dataset
        The remapped dataset in xarray format.
    """
    if remap_method == "bil":
        remap = cdo.Cdo().remapbil(target_grid, input=ds, returnXDataset=True)
    elif remap_method == "con":
        remap = cdo.Cdo().remapcon(target_grid, input=ds, returnXDataset=True)
    elif remap_method == "dis":
        remap = cdo.Cdo().remapdis(target_grid, input=ds, returnXDataset=True)
    elif remap_method == "nn":
        remap = cdo.Cdo().remapnn(target_grid, input=ds, returnXDataset=True)
    else:
        raise ValueError(f"Remap method {remap_method} not supported: choose from 'bil', 'con', 'dis' or 'nn'.")
    if output_path:
        remap.to_netcdf(output_path)
    #Make sure the remap dataset is an dask array
    remap = remap.chunk(chunks="auto")

    #Rename the longitude and latitude to lon and lat if they are not already named like this.
    if "longitude" in remap:
        remap = remap.rename({"longitude": "lon"})
    if "latitude" in remap:
        remap = remap.rename({"latitude": "lat"})
    return remap