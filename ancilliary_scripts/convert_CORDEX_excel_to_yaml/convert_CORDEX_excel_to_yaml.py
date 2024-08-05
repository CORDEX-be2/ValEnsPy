# Convert CORDEX CMIP6 Variables table to yaml for Valenspy

# Inne Vanderkelen, July 2024.

# Original CORDEX file: https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fcordex.org%2Fwp-content%2Fuploads%2F2022%2F09%2FCORDEX_CMIP6_Atmosphere_Variable_List.xlsx&wdOrigin=BROWSELINK

import pandas as pd
import yaml

# identify sheets and header info

sheets = ["Atmos CORE", "Atmos Tier 1", "Atmos Tier 2"]
header = 6

# convert different sheets to ymls
for sheet in sheets:

    # load variable file as a dataframe
    df = pd.read_excel(
        "CORDEX_CMIP6_Atmosphere_Variable_List.xlsx", sheet_name=sheet, header=header
    )
    # create dictionary to write as .yaml
    data_dict = {}
    for _, row in df.iterrows():
        key = row["output variable name"]
        data_dict[key] = {
            "units": row["units"],
            "standard_name": row["standard_name"],
            "long_name": row["long_name"],
            "realm": "atmos",
        }
        if pd.notna(row["ag"]):
            data_dict[key]["aggregation"] = row["ag"]
        if pd.notna(row["Comments"]):
            data_dict[key]["comments"] = row["Comments"]

    # Convert dictionary to YAML and save to file
    with open(sheet + ".yml", "w") as file:
        yaml.dump(data_dict, file, default_flow_style=False)
