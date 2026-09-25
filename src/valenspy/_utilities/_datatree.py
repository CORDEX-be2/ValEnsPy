import pandas as pd
import xarray as xr
from copy import deepcopy

def select_period(dt: xr.DataTree, value: str):
    """All leaves of `dt` that sit under a node named `value`, at any depth - e.g.
    `select_period(dt, "historical")` returns every leaf under any node literally named
    "historical", however deep in the tree. Matches a whole path segment, not a substring
    (a node named "ssp245historical" does NOT match `value="historical"`).

    Unlike `split_by_level`, this needs no fixed tree depth - `dt` can be nested however the
    source ensemble happens to be structured, and different periods don't need to sit at the
    same depth. Used internally by `climate_change_signal_per_member`/`climatology_per_member`/
    `climate_change_signal_ensemble_mean_grid` to select each period by name rather than by a
    project-specific level number.

    Parameters
    ----------
    dt : xr.DataTree
        The DataTree to search.
    value : str
        The node name identifying the period/group to select, e.g. "historical", "ssp245".

    Returns
    -------
    xr.DataTree
        Every leaf of `dt` sitting under a node named `value` - an empty (leafless) DataTree,
        not an error, if no node is named `value` anywhere in `dt`; the caller decides whether
        that's itself an error.
    """
    return dt.filter(lambda node: node.dataset is not None and not node.children and value in node.path.strip("/").split("/"))

def datatree_var_range(dts: xr.DataTree | list, var: str) -> tuple:
    """(min, max) of `var` across every real leaf of `dts` - a single DataTree, or an
    iterable of them (e.g. one per period, combined into one shared range). Leaves without
    `var`, or with no data, are skipped.

    Parameters
    ----------
    dts : xr.DataTree or iterable of xr.DataTree
        Tree(s) to scan.
    var : str
        The variable to compute the range of.

    Returns
    -------
    (float, float)
        (min, max) of `var` across every qualifying leaf.
    """
    if isinstance(dts, xr.DataTree):
        dts = [dts]
    leaves = [leaf.ds[var] for dt in dts for leaf in dt.leaves if leaf.has_data and var in leaf.ds.data_vars]
    if not leaves:
        raise ValueError(f"No leaf has data for variable {var!r}.")
    # Reduce each leaf to a scalar before combining - leaves can sit on different native
    # grids, so concatenating the raw arrays directly is not shape-safe.
    mins = xr.concat([da.min() for da in leaves], dim="_leaf")
    maxs = xr.concat([da.max() for da in leaves], dim="_leaf")
    return float(mins.min()), float(maxs.max())

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