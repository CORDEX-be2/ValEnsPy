# Scripts to Preprocess RADCLIM Dataset

This repository contains scripts to:
1. Convert the RADCLIM data from HDF5 to daily netcdf-files.
2. Submit this process as a job on the VSC HPC.

(c) Kwinten Van Weverberg
August 2024

## Overview of Files

### Scripts
- **`radclim_daily.py`**:
  - Reads in the raw RADCLIM HDF5 files and converts those to daily netcdf-files.
- **`radclim_daily.sh`**:
  - Submits the `radclim_daily.py` as a job on the Tier-1 machine Hortense.

## Usage

0. **Requirements**

    - Tier-1: No additional requirements needed, the necessary modules are automatically loaded by the shell script.
    - kili: The necessary modules can be modified in the shell script to `module load GEOS/3.8.0-GCC-8.3.0-Python-3.7.4`. (untested)
