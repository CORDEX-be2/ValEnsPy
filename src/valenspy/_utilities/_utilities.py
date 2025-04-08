from pathlib import Path
import xarray as xr
from yaml import safe_load
import inspect
import docstring_parser

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

def _set_global_attributes(ds: xr.Dataset, metadata_info):
    for key, value in metadata_info.items():
        ds.attrs[key] = value
    return ds

def _fix_lat_lon(ds: xr.Dataset):
    # rename dimensions if not yet renamed
    if "lon" not in ds.coords:
        ds = ds.rename({"longitude": "lon"})
    if "lat" not in ds.coords:
        ds = ds.rename({"latitude": "lat"})

    # make sure lat and lon are sorted ascending
    ds = ds.sortby("lat").sortby("lon")
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
    #If yml_name is a path
    if isinstance(yml_name, Path):
        file = yml_name
    else:
        src_path = Path(__file__).resolve().parent.parent
        file = src_path / "ancilliary_data" / f"{yml_name}.yml"
    with open(file, "r") as file:
        yml = safe_load(file)
    return yml

def generate_parameters_doc(func):
        """Generate the parameters section of the docstring to match the signature of the diagnostic function."""

        signature = inspect.signature(func)
        docstring = inspect.getdoc(func)
        doc_lines = ["Parameters", "----------"]

        parsed = docstring_parser.parse(docstring) if docstring else None
        param_docs = {p.arg_name: p.description for p in parsed.params} if parsed else {}

        for name, param in signature.parameters.items():
            # Handle special parameter types (*args, **kwargs)
            display_name = name
            if param.kind == inspect.Parameter.VAR_POSITIONAL:
                display_name = f"*{name}"
            elif param.kind == inspect.Parameter.VAR_KEYWORD:
                display_name = f"**{name}"

            # Get annotation
            if param.annotation is inspect.Parameter.empty:
                annotation_str = "Any"
            elif hasattr(param.annotation, "__name__"):
                annotation_str = param.annotation.__name__
            else:
                annotation_str = str(param.annotation)

            # Get default
            if param.default is not inspect.Parameter.empty:
                default_str = f", default={param.default!r}"
            else:
                default_str = ""

            # Description (from docstring or fallback)
            description = param_docs.get(name, f"Description of {name}.")

            doc_lines.append(f"{display_name} : {annotation_str}{default_str}")
            doc_lines.append(f"    {description}")

        return "\n".join(doc_lines) + "\n\n"