from scipy.spatial import cKDTree
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable

def point_data_2_climate_grid(df_obs, ds_mod, variable):
    """Averages observational data sets of which some observations are located within one grid cell of the model domain.

    Parameters
    ----------
    df_obs : dataframe that contains the observational data
    ds_mod: dataset that contains the model data.
    variable: list of variables that should be averaged

    Returns
    -------
    pd.DataFrame
        DataFrame in which the observations that lie within the same model gridcell are averaged.
    """

    if isinstance(variable, str):
        variable = [variable]

    # Get lat/lon from the model dataset (assuming it's structured with 2D lat/lon variables)
    model_lats = ds_mod['lat'].values  # Extract latitudes
    model_lons = ds_mod['lon'].values  # Extract longitudes

    # Flatten the lat/lon arrays in case they are 2D
    model_coords = np.column_stack([model_lats.ravel(), model_lons.ravel()])

    # Build a KDTree for fast nearest-neighbor lookup
    tree = cKDTree(model_coords)

    # Query nearest grid point for each observation
    obs_coords = df_obs[['lat', 'lon']].values
    distances, indices = tree.query(obs_coords)
    
    # Assign observations to grid cell indices
    df_obs['grid_index'] = indices
    df_indices= df_obs[["lat", "lon", "grid_index"]].drop_duplicates().reset_index(drop = True)

    # Compute mean values for each grid cell
    df_columns = variable
    df_columns.extend(["grid_index", "time"])
    obs_aggregated = df_obs[df_columns].groupby(["time", "grid_index"]).mean()

    # Merge with grid coordinates for reference
    df_obs['grid_index'] = indices
    df_indices= df_obs[["lat", "lon", "grid_index"]].drop_duplicates().reset_index(drop = True)
    df_indices = df_indices.rename({"orig_lat" : "lat", "orig_lon": "lon"})

    # Reset index to expose "grid_index" as a column
    obs_aggregated = obs_aggregated.reset_index()

    # Now correctly assign lat/lon from model_coords using grid_index
    obs_aggregated["lat"] = model_coords[obs_aggregated["grid_index"], 0]
    obs_aggregated["lon"] = model_coords[obs_aggregated["grid_index"], 1]

    return [obs_aggregated, df_indices]


def plot_points_map(d_coord_points, region, bounds, scale_var=None, lat_name = "lat", lon_name = "lon", station_id = "code", add_labels = False): 
    """
    Plots a map with points and an optional colorbar.
    
    Parameters:
    - d_coord_points (pd.DataFrame): DataFrame containing 'lon', 'lat', 'code', and optionally scale_var.
    - region (str): The region to display.
    - bounds (dict): Dictionary containing lat/lon bounds for regions.
    - scale_var (str, optional): Column name for variable used to scale colors.
    """
    # Create a figure and set the projection to PlateCarree
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Add coastlines and country borders
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # Handle scaling variable
    if scale_var:
        cmap = plt.cm.viridis
        norm = Normalize(vmin=d_coord_points[scale_var].min(), vmax=d_coord_points[scale_var].max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)  # Create a ScalarMappable for the colorbar
        sm.set_array([])  # ScalarMappable requires a linked array but it won't be used here

    # Plot each point and add a label
    for i in np.arange(0, len(d_coord_points)):
        lon = d_coord_points.loc[i, lon_name]
        lat = d_coord_points.loc[i, lat_name]
        key = d_coord_points.loc[i, station_id]

        if scale_var:
            var_value = d_coord_points.loc[i, scale_var]
            code_color = cmap(norm(var_value))
            ax.scatter(lon, lat, edgecolor='black', s=50, c=[code_color], transform=ccrs.PlateCarree())
        else:
            ax.scatter(lon, lat, color='red', s=50, transform=ccrs.PlateCarree())
        
        # Add text labels
        if add_labels:
            ax.text(lon + 0.05, lat - 0.05, key, transform=ccrs.PlateCarree())


    # Set extent based on bounds
    lat_bounds = bounds[region]['lat_bounds']
    lon_bounds = bounds[region]['lon_bounds']
    ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())

    # Add colorbar if scale_var is provided
    if scale_var:
        cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.02, shrink=0.7)
        cbar.set_label(f'{scale_var.capitalize()}')

    # Add a title
    ax.set_title('Location of Points', loc='right')

    # Show the plot
    plt.show()

    return ax