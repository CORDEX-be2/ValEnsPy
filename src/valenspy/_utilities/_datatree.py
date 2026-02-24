import pandas as pd
import xarray as xr
from copy import deepcopy

def split_by_level(dt: xr.DataTree, level: int):
    """
    Split a DataTree into multiple DataTrees based on the unique values at a given level in the node paths.

    Parameters
    ----------
    dt : xr.DataTree
        The DataTree to split.
    level : int
        The level in the node paths to split the DataTree by. Level 0 is the root level.
    
    Returns
    -------
    dict
        A dictionary where the keys are the unique values at the specified level in the node paths, and the values are the corresponding DataTrees containing only the nodes with that value at the specified level.
    """
    #Make a deep copy of the datatree to avoid modifying the original one. This is needed as we are going to orphan the datasets in the new datatrees, which would also orphan them in the original datatree if we don't make a copy.
    dt = deepcopy(dt)
    if level >= 2:
        dt = restructure_by_level(dt, level)
    result = {value: dt[value] for value in set(dt.children)}
    for value in result:
        result[value].orphan()
    return result

def restructure_by_level(dt: xr.DataTree, level: int):
    """
    Restructure a DataTree such that level n in the node paths becomes the new root level.

    Parameters
    ----------
    dt : xr.DataTree
        The DataTree to restructure.
    level : int
        The level in the node paths to restructure the DataTree by. Level 0 is the root level.

    Returns
    -------
    xr.DataTree
        A restructured DataTree where level n in the node paths becomes the new root level.
    """
    if level < 2:
        raise ValueError("Level must be greater than or equal to 2 as reording the first level leaves the tree unchanged")
    level = level - 1 #As the root level is not included in the path split, we need to subtract 1 from the level to get the correct index.
    reorganized_nodes = {
        "/".join([path.split("/")[level]] + path.split("/")[:level] + path.split("/")[level+1:]): node.dataset
        for path, node in dt.subtree_with_keys
        if len(path.split("/")) > level #Strict
    }
    return xr.DataTree.from_dict(reorganized_nodes)

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