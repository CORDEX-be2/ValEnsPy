import pandas as pd
import xarray as xr
from copy import deepcopy

def mask_to_reference_coverage(dt: xr.DataTree, reference: xr.Dataset, dim: str = "time") -> xr.DataTree:
    """Mask every leaf of `dt` to `reference`'s own coverage: NaN wherever `reference`
    has NO valid value anywhere along `dim`, at every OTHER coordinate. Useful e.g. for
    restricting model data to the same spatial extent an observational reference actually
    covers (so a comparison isn't biased by grid cells - like ocean pixels for a
    land-only observational dataset - the reference never has real data for), with
    `dim="time"` (the default) checking "ever valid at any timestep" rather than
    requiring validity at every timestep.

    Parameters
    ----------
    dt : xr.DataTree
        The data to mask - every leaf is restricted independently, against the SAME
        `reference` coverage.
    reference : xr.Dataset
        The reference dataset defining "valid coverage" - `reference.notnull().any(dim=
        dim)` gives the mask (True wherever `reference` has at least one valid value
        along `dim`), applied to every leaf of `dt` via `.where(...)`. Typically already
        restricted to the same variable(s) as `dt`, so coverage reflects only what's
        actually being compared.
    dim : str, optional
        The dimension to reduce over when checking coverage - "does `reference` have a
        valid value anywhere along this dimension". Default "time" (mask by spatial
        extent, checking temporal validity); pass e.g. `dim=["lat", "lon"]` instead to
        mask by TEMPORAL coverage (checking spatial validity) if that's what's needed.

    Returns
    -------
    xr.DataTree
        `dt`, masked to `reference`'s own coverage.
    """
    coverage = reference.notnull().any(dim=dim)

    def _mask(ds):
        if not ds:  # Empty leaf placeholder - nothing to mask.
            return ds
        return ds.where(coverage)

    return dt.map_over_datasets(_mask)

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

def reorder(dt: xr.DataTree, new_order: list):
    """
    Reorder the datatree paths according to the order specified in new_order which should be a list of all leave paths in the datatree. 
    The index of the leave order will be used as a new parent node to order the datatree by. 

    Parameters
    ----------
    dt : xr.DataTree
        The DataTree to reorder.
    new_order : list
        A list of all leave paths in the datatree in the desired order. The index
        of the leave order will be used as a new parent node to order the datatree by.
    
    Returns
    -------
    xr.DataTree
        A reordered DataTree where the paths are reordered according to the order specified in new_order.
    """
    return xr.DataTree.from_dict({
        f"{i}/{path}" : dt[path].dataset
        for i, path in enumerate(new_order)
    })

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

def restructure_by_attributes(dt: xr.DataTree, attributes: list):
    """
    Restructure a DataTree such that the specified attributes form the new path of the DataTree. The attributes should be specified in the order they should appear in the new path.

    Parameters
    ----------
    dt : xr.DataTree
        The DataTree to restructure.
    attributes : list
        A list of attributes to restructure the DataTree by. The attributes should be specified in the order they should appear in the new path.

    Returns
    -------
    xr.DataTree
        A restructured DataTree where the specified attributes form the new path of the DataTree.
    """
    reorganized_nodes = {
        "/".join([str(node.dataset.attrs.get(attr, "None")) for attr in attributes]): node.dataset
        for path, node in dt.subtree_with_keys
        if path and (node.dataset is not None) and (not node.children) #Only include leaf nodes with datasets
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
                attr_dict = {}
                if isinstance(add_attributes, bool): #If add_attributes is True, we add all attributes. If it's a list, we only add the specified attributes.
                    attr_dict.update(ds.attrs) #The dataset attributes
                    for var in ds.data_vars: #The variable attributes
                        attr_dict.update({f"{var}_{attr}": ds[var].attrs[attr] for attr in ds[var].attrs})
                else:
                    for attr in add_attributes:
                        #if attr is a substring of any attribute in ds.attrs:
                        for ds_attr in ds.attrs:
                            if attr in ds_attr:
                                attr_dict[attr] = ds.attrs[ds_attr]
                                break #Break the loop after finding the first match to avoid adding multiple attributes with the same substring.
                #If the attribute is not a string try to cast it to a string, else drop it
                for attr in attr_dict:
                    if not isinstance(attr_dict[attr], str):
                        try:
                            attr_dict[attr] = str(attr_dict[attr])
                        except:
                            del attr_dict[attr]
                df = df.assign(**attr_dict)
            data_frames.append(df)
                
    return pd.concat(data_frames, axis=0).reset_index()