# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
if __name__ == '__main__':
    import valenspy
    import seaborn as sns
    import xarray as xr
    import matplotlib.pyplot as plt
    from xclim.indicators.atmos import cooling_degree_days

    variables = ["pr","tas", "tasmax", "tasmin"]
    period=["1980","2020"]
    frequency=["day","daily"]

    time_decorder = xr.coders.CFDatetimeCoder() #Dealing with different calendars

    xarray_open_kwargs = {
    "chunks": {"time": 1},
    "decode_times": time_decorder,
    "decode_coords": "all",
    }

    # %%
    #Creating the input manager and adding our personal yml file to the catalog
    m = valenspy.InputManager("bilan") #This takes very long for the all the CORDEX data...  

    # %%
    cat = m.search(
    source_type=["CMIP6"],
    experiment_id=["historical"],
    time_period=["2000","2001"]
    )
    cat.unique()

    # %%
    # Set the scheduler to 'processes' for parallel processing.
    import dask
    dask.config.set(scheduler='processes')

    import os
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    ids = m.unique()['source_id']
    ids.remove("CESM2-WACCM-FV2") #Removed froom the list due to not having monotonic global attributes

    dt = m.search(
    source_type=["CMIP6"],
    experiment_id=["historical"],
    source_id=ids,
    time_period=period,
    variable_id=["pr"],
    ).to_datatree(
    xarray_open_kwargs=xarray_open_kwargs,
    xarray_combine_by_coords_kwargs={
        "coords": "different",
    },
    levels=["source_id", "realization"]
    )

    # %%
    import xclim
    dt_pr_days = valenspy.xclim_indicator(dt, xclim.indices.max_1day_precipitation_amount, vars=["pr"])
    dt_pr_days = valenspy.convert_units_to(dt_pr_days, var="pr", target_unit="mm/day")

    Ukkel = (4.37, 50.79)
    dt_pr_days = dt_pr_days.sel(lat=Ukkel[1], lon=Ukkel[0], method="nearest")

    # %%
    from dask.diagnostics import ProgressBar
    with ProgressBar():
        dt_pr_days = dt_pr_days.compute()

    # %%
    import pandas as pd
    data_frames = []
    for key, ds in dt_pr_days.to_dict().items():
        if ds:
            df = ds.to_dataframe()
            df["source_id"] = key
            data_frames.append(df)
    df = pd.concat(data_frames)
    df = df.reset_index()
    # df["time"] = pd.to_datetime(df["time"])


    # %%
    plt.figure(figsize=(12, 6))
    sns.boxplot(y="pr", data=df, hue="source_id", palette="Set2")
    plt.title("Maximum 1-day Precipitation Amount (mm/day)")
    #Save the figure
    plt.savefig("max_1day_precipitation_amount.png", dpi=300, bbox_inches="tight")

    # %%
    # df["time"] = df["time"].dt.year
    #Convert to the same calendar format

    # sns.barplot(data=df,x="time",y="pr", hue="source_id", palette="Set2")
    # ## Make the x axis only the year 

    # plt.xticks(rotation=90)
    # plt.ylabel("Maximum 1-day Precipitation Amount (mm/day)")

    # plt.xlabel("Year")
