"""
Author: K. Van Weverberg, August 2024
"""

import h5py
import netCDF4 as nc
import numpy as np
import datetime
import os
from pyproj import Proj, transform
import zipfile


def estimate_corners(lat, lon):
    """
    Estimate the corners of each grid cell for 2D curvilinear grids.
    The corners are estimated by averaging the midpoints of neighboring cells.
    Args:
        lat (np.array): 2D array of latitudes for cell midpoints.
        lon (np.array): 2D array of longitudes for cell midpoints.
    Returns:
        lat_corners, lon_corners: 2D arrays of shape (nlat, nlon, 4) containing the estimated
                                  latitude and longitude of each corner for each grid cell.
    """
    nlat, nlon = lat.shape
    lat_corners = np.zeros((nlat, nlon, 4))  # Initialize arrays for corners
    lon_corners = np.zeros((nlat, nlon, 4))

    # Internal cells
    for i in range(0, nlat):
        for j in range(0, nlon):
            im1 = np.max([0, i - 1])
            ip1 = np.min([nlat - 1, i + 1])
            jm1 = np.max([0, j - 1])
            jp1 = np.min([nlon - 1, j + 1])

            # Average to get the SW corner of the (i, j) cell
            lat_corners[i, j, 0] = (
                lat[i, j] + lat[i, jm1] + lat[im1, j] + lat[im1, jm1]
            ) / 4
            lon_corners[i, j, 0] = (
                lon[i, j] + lon[i, jm1] + lon[im1, j] + lon[im1, jm1]
            ) / 4

            # Average to get the NW corner of the (i, j) cell
            lat_corners[i, j, 1] = (
                lat[i, j] + lat[i, jp1] + lat[im1, j] + lat[im1, jp1]
            ) / 4
            lon_corners[i, j, 1] = (
                lon[i, j] + lon[i, jp1] + lon[im1, j] + lon[im1, jp1]
            ) / 4

            # Average to get the NE corner of the (i, j) cell
            lat_corners[i, j, 2] = (
                lat[i, j] + lat[i, jp1] + lat[ip1, j] + lat[ip1, jp1]
            ) / 4
            lon_corners[i, j, 2] = (
                lon[i, j] + lon[i, jp1] + lon[ip1, j] + lon[ip1, jp1]
            ) / 4

            # Average to get the SE corner of the (i, j) cell
            lat_corners[i, j, 3] = (
                lat[i, j] + lat[i, jm1] + lat[ip1, j] + lat[ip1, jm1]
            ) / 4
            lon_corners[i, j, 3] = (
                lon[i, j] + lon[i, jm1] + lon[ip1, j] + lon[ip1, jm1]
            ) / 4

    return lat_corners, lon_corners


def read_hdf5_radar_data(hdf5_file_path):
    """
    Read radar data from an HDF5 file and extract relevant information
    following the RADQPE User Guide specifications.
    """
    with h5py.File(hdf5_file_path, "r") as hdf:
        # Assuming dataset1 contains the data of interest
        dataset_path = "/dataset1/data1/data"
        data = hdf[dataset_path][:]

        # Extracting geolocation and time information
        geoloc_attrs = hdf["/dataset1/where"].attrs
        geoloc_attrs_proj = hdf["/where"].attrs
        time_attrs = hdf["/dataset1/what"].attrs

        ll_lat = hdf["where"].attrs["LL_lat"]
        ll_lon = hdf["where"].attrs["LL_lon"]
        lr_lat = hdf["where"].attrs["LR_lat"]
        lr_lon = hdf["where"].attrs["LR_lon"]
        ul_lat = hdf["where"].attrs["UL_lat"]
        ul_lon = hdf["where"].attrs["UL_lon"]
        ur_lat = hdf["where"].attrs["UR_lat"]
        ur_lon = hdf["where"].attrs["UR_lon"]
        c_lat = hdf["where"].attrs["lat"]
        c_lon = hdf["where"].attrs["lon"]
        projdef = hdf["where"].attrs["projdef"]

        ul_x = hdf["dataset1"]["where"].attrs["UL_x"]
        ul_y = hdf["dataset1"]["where"].attrs["UL_y"]
        xsize = hdf["dataset1"]["where"].attrs["xsize"]
        ysize = hdf["dataset1"]["where"].attrs["ysize"]
        xscale = hdf["dataset1"]["where"].attrs["xscale"]
        yscale = hdf["dataset1"]["where"].attrs["yscale"]

        lr_x = ul_x + (xsize * xscale)
        lr_y = ul_y - (ysize * yscale)

        x = np.arange(ul_x, lr_x, xscale) + xscale / 2
        y = np.arange(lr_y, ul_y, yscale) - yscale / 2

        xleft = np.arange(ul_x, lr_x, xscale)
        xright = np.arange(ul_x, lr_x, xscale) + xscale

        ytop = np.arange(lr_y, ul_y, yscale)
        ybot = np.arange(lr_y, ul_y, yscale) - yscale

        xx, yy = np.meshgrid(x, y)

        yy = np.flip(yy)

        inProj = Proj(
            r"+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
        )
        outProj = Proj(init="epsg:4326")

        lonrad, latrad = transform(inProj, outProj, xx, yy)

        time_info = {
            "startdate": time_attrs["startdate"].decode("utf-8"),
            "starttime": time_attrs["starttime"].decode("utf-8"),
            "enddate": time_attrs["enddate"].decode("utf-8"),
            "endtime": time_attrs["endtime"].decode("utf-8"),
        }

        return data, lonrad, latrad, time_info


start_dt = "2017-01-01 01:00:00"
end_dt = "2022-12-31 23:00:00"

output_directory = "/data/gent/vo/002/gvo00202/rcmgeo/vsc45302/"  # Update this path
directory = "/data/gent/vo/002/gvo00202/rcmgeo/vsc45302/"  # Update this path

dt_format = "%Y%m%d%H%M%S"
start_timestamp = datetime.datetime.strptime(start_dt, "%Y-%m-%d %H:%M:%S")
end_timestamp = datetime.datetime.strptime(end_dt, "%Y-%m-%d %H:%M:%S")
current_timestamp = start_timestamp

data_list = []
time_info_list = []

# Generates file paths for radar data files within a specified date/time range.
while current_timestamp <= end_timestamp:
    year = current_timestamp.strftime("%Y")
    month = current_timestamp.strftime("%m")
    day = current_timestamp.strftime("%d")
    hour = current_timestamp.strftime("%H")
    print(hour)

    filenamezip = os.path.join(
        directory,
        "mnt/HDS_RADCLIM/RADCLIM/export",
        year,
        month,
        day,
        "%s%s%s.radclim.1h.hdf.zip" % (year, month, day),
    )
    with zipfile.ZipFile(filenamezip, "r") as zip_ref:
        #       zip_ref.extractall(directory)

        for member in zip_ref.namelist():
            # Extract the file name (ignoring directories)
            filename = os.path.basename(member)
            if not filename:
                continue  # Skip directories

            # Define the full path for the extracted file
            source = zip_ref.open(member)
            target_path = os.path.join(
                directory, "mnt/HDS_RADCLIM/RADCLIM/export", year, month, day, filename
            )

            # Write the content of the file
            with open(target_path, "wb") as target:
                target.write(source.read())

    # zipfile.ZipFile(filenamezip, 'r')
    filename = current_timestamp.strftime(dt_format) + ".radclim.1h.hdf"
    file_path = os.path.join(
        directory, "mnt/HDS_RADCLIM/RADCLIM/export", year, month, day, filename
    )

    data, lons, lats, time_info = read_hdf5_radar_data(file_path)

    data_list.append(data)
    # Extract time info from filename
    base_filename = os.path.basename(filename)
    file_timestamp = base_filename[:14]  # Assuming yyyymmddhhmmss format
    time_info_list.append(file_timestamp)

    current_timestamp += datetime.timedelta(hours=1)
    next_timestamp = current_timestamp + datetime.timedelta(hours=1)
    nextmonth = next_timestamp.strftime("%m")

    if hour == "23":
        data_array = np.stack(data_list, axis=0)  # Stack along new axis (time)

        nx = data_array.shape[2]
        ny = data_array.shape[1]

        lat_bounds, lon_bounds = estimate_corners(lats, lons)

        netcdf_file_path = os.path.join(
            output_directory, "%s%s%s.radclim.1h.nc" % (year, month, day)
        )

        # Creating a NetCDF file
        ncfile = nc.Dataset(netcdf_file_path, "w", format="NETCDF4")

        # Create dimensions
        lat_dim = ncfile.createDimension("nlat", data_array.shape[1])
        lon_dim = ncfile.createDimension("nlon", data_array.shape[2])
        ncfile.createDimension("nv", 4)
        time = ncfile.createDimension("time", data_array.shape[0])

        # Create variables
        lat_var = ncfile.createVariable("lat", np.float32, ("nlat", "nlon"))
        lon_var = ncfile.createVariable("lon", np.float32, ("nlat", "nlon"))
        lat_bounds_var = ncfile.createVariable(
            "lat_bounds", np.float32, ("nlat", "nlon", "nv")
        )
        lon_bounds_var = ncfile.createVariable(
            "lon_bounds", np.float32, ("nlat", "nlon", "nv")
        )
        times = ncfile.createVariable("time", "f8", ("time",))
        precipitation = ncfile.createVariable(
            "precipitation",
            "f4",
            (
                "time",
                "nlat",
                "nlon",
            ),
            fill_value=-9999.0,
        )

        lat_bounds_var[:, :, :] = lat_bounds
        lon_bounds_var[:, :, :] = lon_bounds

        # Now add the bounds attribute to the 'lat' and 'lon' variables
        ncfile.variables["lat"].bounds = "lat_bounds"
        ncfile.variables["lon"].bounds = "lon_bounds"

        # Optionally, add attributes to the boundary variables
        lat_bounds_var.units = "degrees_north"
        lon_bounds_var.units = "degrees_east"

        # Assign data
        lat_var[:] = lats
        lon_var[:] = lons
        precipitation[:, :, :] = data_array
        # Assigning time data
        dt_format = "%Y%m%d%H%M%S"
        times_list = [
            datetime.datetime.strptime(ts, dt_format) for ts in time_info_list
        ]
        times.units = "hours since 2017-01-01 00:00:00"
        times.calendar = "gregorian"
        times[:] = nc.date2num(times_list, units=times.units, calendar=times.calendar)

        # Add attributes for CF compliance
        lat_var.units = "degrees_north"
        lat_var.standard_name = "latitude"
        lon_var.units = "degrees_east"
        lon_var.standard_name = "longitude"
        precipitation.units = (
            "mm"  # Assuming radar data is precipitation; adjust as needed
        )
        precipitation.long_name = "Radar estimated precipitation"
        precipitation.coordinates = "lat lon"

        # Global attributes
        ncfile.title = "Radar Data with Curvilinear Grid"
        ncfile.source = "Generated for demonstration"
        ncfile.Conventions = "CF-1.8"

        ncfile.close()

        data_list = []
        time_info_list = []
