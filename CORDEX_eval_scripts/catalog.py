"""
catalog.py

Adapted from the original script by the EURO-CORDEX joint evaluation project:

Functions:
- human_readable(df): Creates a human-readable summary of the dataset.
- create_excel(filename): Creates a human-readable Excel file from the dataset.
- update_catalog(catalog, root): Updates the catalog with metadata from the specified root directory.
"""

import os
import re
import pandas as pd
from os import path as op

root_dic = {
    "ALARO": "/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/CORDEXbeII/run_ALARO_sfx/out/PROD_ERA5_EUR11/19780101/CORDEX/CMIP6/",
}

CATALOG = "catalog.csv"

COLS = [
    "project_id",
    "mip_era",
    "activity_id",
    "domain_id",
    "institution_id",
    "driving_source_id",
    "driving_experiment_id",
    "driving_variant_label",
    "source_id",
    "version_realization",
    "frequency",
    "version",
    "time_range",
    "variable_id",
]


def create_path_pattern(drs, sep="/"):
    attrs = drs.split(sep)
    drs = sep.join([f"(?P{attr}[^/]+)" for attr in attrs])
    # Allow for an optional root directory
    drs = r"^/?(?:[^/]+/)*" + drs
    return re.compile(drs)


def parse_filepath(filename, path_format):
    # pattern = create_pattern(drs)
    if path_format == "cmip6-cordex":
        regex = r"(?P<project_id>[^/]+)/(?P<mip_era>[^/]+)/(?P<activity_id>[^/]+)/(?P<domain_id>[^/]+)/(?P<institution_id>[^/]+)/(?P<driving_source_id>[^/]+)/(?P<driving_experiment_id>[^/]+)/(?P<driving_variant_label>[^/]+)/(?P<source_id>[^/]+)/(?P<version_realization>[^/]+)/(?P<frequency>[^/]+)/(?P<variable_id>[^/]+)/(?P<version>[^/]+)/(?P<filename>(?P<variable_id_2>[^_]+)_(?P<domain_id_2>[^_]+)_(?P<driving_source_id_2>[^_]+)_(?P<driving_experiment_id_2>[^_]+)_(?P<driving_variant_label_2>[^_]+)_(?P<institution_id_2>[^_]+)_(?P<source_id_2>[^_]+)_(?P<version_realization_2>[^_]+)_(?P<frequency_2>[^_]+)(?:_(?P<time_range>[^.]+))?\.nc)"
        regex = r"^/?(?:[^/]+/)*" + regex
    pattern = re.compile(regex)
    match = pattern.match(filename)
    if match:
        return match.groupdict()
    else:
        raise ValueError("The filepath does not match the expected pattern.")

def parse_catalog_files(root, path_format):
    datasets = []
    # Define the regex pattern for the filename
    for root, dirs, files in os.walk(root):
        # only parse if files found
        if not files:
            continue
        for file in files:
            if ".nc" in file:
                filename = op.join(root, file)
                print(f"parsing {filename}")
                metadata = parse_filepath(filename, path_format)
                metadata["path"] = filename
                datasets.append(metadata)
    return datasets


def human_readable(df):
    """
    Creates a human-readable summary of the dataset.

    Parameters:
    df (pandas.DataFrame): The input DataFrame containing the dataset.

    Returns:
    pandas.DataFrame: A DataFrame with grouped and summarized data.
    """
    cols = [item for item in COLS if item not in ["variable_id", "time_range"]]

    def to_list(x):
        return list(dict.fromkeys(list(x)))

    return df.groupby(cols)["variable_id"].apply(to_list).to_frame()  # .reset_index()


def create_catalog(root, path_format):
    """
    Updates the catalog with metadata from the specified root directory.

    Parameters:
    root (str): The root directory to scan for metadata.
    path_format (str): The format of the paths to parse.

    Returns:
    pandas.DataFrame: A DataFrame containing the most recent catalog.
    """
    df = pd.DataFrame(parse_catalog_files(root, path_format))[COLS + ["path"]]
    return df


if __name__ == "__main__":
    df_ALARO = create_catalog(root_dic["ALARO"], "cmip6-cordex")
    df = pd.concat([df_ALARO])
    df.to_csv(f"{CATALOG}", index=False)

    human_readable(df).to_csv(f"{CATALOG.replace('.csv', '_summary.csv')}")