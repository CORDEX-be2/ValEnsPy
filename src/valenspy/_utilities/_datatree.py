import pandas as pd
import xarray as xr

def datatree_to_dataset(dt: xr.DataTree, **kwargs):
    """
    Convert a DataTree to a xarray Dataset.

    Parameters
    ----------
    dt : xr.DataTree
        The DataTree to convert to a xarray Dataset.
    **kwargs : dict
        Keyword arguments to pass to the xarray concat function.
    """
    datasets = []
    for key, ds in dt.to_dict().items():
        if ds:
            ds_copy = ds.copy()
            ds_copy = ds_copy.expand_dims({"id": [str(key)]})
            datasets.append(ds_copy)
            
    return xr.concat(datasets, dim="id", **kwargs)
    
def datatree_to_dataframe(dt: xr.DataTree, add_attributes=False):
    """
    Convert a DataTree to a pandas DataFrame.
    """
    data_frames = []
    for key, ds in dt.to_dict().items():
        if ds:
            if not ds.dims:  # Non-dimensional datasets
                df = pd.DataFrame({var: float(ds[var].values) for var in ds.data_vars}, index=[key])
            else:  # Dimensional datasets
                df = ds.to_dataframe()
            df["id"] = str(key)
            
            if add_attributes:
                for attr in add_attributes:
                    #if attr is a substring of any attribute in ds.attrs:
                    for ds_attr in ds.attrs:
                        if attr in ds_attr:
                            df[attr] = ds.attrs[ds_attr]
                            break
            data_frames.append(df)
                
    return pd.concat(data_frames, axis=0).reset_index()