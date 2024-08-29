from pathlib import Path
import xarray as xr
from yaml import safe_load


def load_xarray_from_data_sources(data_sources):
    """Return a xarray dataset from an list of input files (Path objects), a single input file (Path object) or a xarray dataset.
    This utility function enables the user to input different types of inputs which are then converted to a xarray dataset.

    Parameters
    ----------
    inputs : Path or list(Path) or xarray.Dataset
        The input file or list of input files to convert.

    Returns
    -------
    xarray.Dataset
        An xarray dataset.
    """

    if isinstance(data_sources, Path) or isinstance(data_sources, list):
        ds = xr.open_mfdataset(data_sources, combine="by_coords", chunks="auto")
    elif isinstance(data_sources, xr.Dataset):
        ds = data_sources
    else:
        raise TypeError(
            "The input should be a Path, a list of Paths or an xarray dataset."
        )
    return ds


def load_yml(yml_name):
    """Load a yaml file into a dictionary from the ancilliary_data folder. The yaml file should be in the ancilliary_data folder.

    Parameters
    ----------
    yml_name : str
        The name of the yaml file to load.

    Returns
    -------
    dict
        The yaml file loaded into a dictionary.
    """
    src_path = Path(__file__).resolve().parent.parent
    with open(src_path / "ancilliary_data" / f"{yml_name}.yml") as file:
        yml = safe_load(file)
    return yml
