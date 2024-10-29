# Scripts to download and preprocess ERA5 and ERA5-Land

This repository contains a collection of Python scripts to:
1. Download the raw hourly ERA5 and ERA5-Land dataset from the CDS
2. Merge the datasets into daily and monthly timesteps

(c) Inne Vanderkelen, Bert van Schaeybroeck
August 2024

## Overview of Files

### Scripts
- **`download_era5_single_level.py`**:
  - Connects to the CDS and downloads the ERA5 data. Variables and time period to be determined by the user. 
- **`process_era5_hourly_to_daily_monthly.py`**:
  - Merge the datasets into daily and monthly timesteps if not yet available

## Usage

0. **Requirements**

    CDO module loaded and the valenspy dev python environment and additionaly, the oracledb python package
  
1. **Download and Preprocess the Data:**
   ```bash
   python process_era5_hourly_to_daily_monthly
    ```
