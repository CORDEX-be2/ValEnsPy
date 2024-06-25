# Scripts to Preprocess CLIMATE_GRID Dataset

This repository contains a collection of Python scripts to:
1. Download the CLIMATE_GRID dataset from the RMI Oracle database as vector data and store it as CSV.
2. Convert the dataset to NetCDF format with the necessary attributes.
3. Remap the dataset to a self-defined lat-lon grid.

(c) Inne Vanderkelen, Bert van Schaeybroeck, Nicolas Ghilain  
June 2024

## Overview of Files

### Scripts
- **`preprocess_CLIMATE_GRID.py`**: 
  - Connects to the RMI Oracle database and downloads the CLIMATE_GRID data.
  - Converts and saves the data as a NetCDF file.
- **`remap_CLIMATE_GRID_to_latlon.py`**: 
  - Remaps the projected (raw) CLIMATE_GRID NetCDF data to a lat/lon grid.
  - The target grid can be user-defined within the script.
- **`preprocess_CLIMATE_GRID_Rdata_to_netcdf.py`**: 
  - Converts the Rdata downloaded from the Oracle database using Bert's script to a NetCDF file.

### Data Files
- **`grid_5kmx5km.csv`**: 
  - Contains pixel lat and lon variables to transpose to a user-defined grid.
  - Used by `preprocess_CLIMATE_GRID.py`.
- **`CLIMATE_GRID_meta.csv`**: 
  - Contains the variable names and units for the CLIMATE_GRID dataset.
  - Used to write the NetCDF file by `preprocess_CLIMATE_GRID.py`.
- **`lambert_coordinates_full_climate_grid.csv`**: 
  - Lists pixel IDs and corresponding x, y coordinates in LAMBERT2008.
  - Used by `preprocess_CLIMATE_GRID.py`.

## Usage

0. **Requirements**

    CDO module loaded and the valenspy dev python environment and additionaly, the oracledb python package 
    ```bash 
    conda install oracledb 
    ```

1. **Download and Preprocess the Data:**
   ```bash
   python preprocess_CLIMATE_GRID.py
    ```