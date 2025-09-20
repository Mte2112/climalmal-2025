import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
import sys


def fliplon(ds):
    """
    flip lon from 0-360 to -180-180
    """

    if "lon" in ds.dims:
        lon = "lon"
        lat = "lat"
    else:
        lon = "longitude"
        lat = "latitude"

    # make sure lon/lat indices are formatted in order for slicing (ndvi lat is backwards)
    ds = ds.sortby(["time", lat, lon])

    # adjust lon values to make sure they are within (-180, 180)
    # allows for easier plotting and compatibility with certain packages
    ds["_longitude_adjusted"] = xr.where(ds[lon] > 180, ds[lon] - 360, ds[lon])

    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (
        ds.swap_dims({lon: "_longitude_adjusted"})
        .sel(**{"_longitude_adjusted": sorted(ds._longitude_adjusted)})
        .drop(lon)
    )

    ds = ds.rename({"_longitude_adjusted": lon})

    return ds


def process_nc_malawi(path, year_start, year_end, mal_bounds, operation):
    """
    load and reformat netcdf data, slice malawi bounds
    """

    ds = xr.open_mfdataset(path)
    # make sure lon/lat indices are formatted in order for slicing (ndvi lat is backwards)
    ds = ds.sortby(["time", "lat", "lon"])

    # flip lon -180 to 180
    ds = fliplon(ds)

    # get box around malawi
    ds_mal = ds.sel(
        lon=slice(mal_bounds[0], mal_bounds[1]), lat=slice(mal_bounds[2], mal_bounds[3])
    )
    ds_mal = ds_mal.sel(time=slice(str(year_start), str(year_end)))
    ds_mal = ds_mal.sortby("time")
    if operation == "mean":
        ds_mal = ds_mal.resample(time="1M").mean()  # get monthly mean
    elif operation == "sum":
        ds_mal = ds_mal.resample(time="1M").sum()  # get monthly sum
    elif operation == "raw":
        None  # keep data in native format
    else:
        print("operation not recognized. exiting...")
        sys.exit()
    # change all the days to 15th to match the incidence data format
    time_index = pd.to_datetime(ds_mal["time"].values)
    new_times = pd.to_datetime([f"{t.year}-{t.month}-15" for t in time_index])
    ds_mal["time"] = xr.DataArray(new_times.values, dims="time")

    return ds_mal


def process_malaria_data(
    path, year_start, year_end, age_group, dist_to_remove, columns_to_keep
):
    """
    process malaria incidence data
    """

    # read in file
    mal = pd.read_csv(path)
    if columns_to_keep:
        mal = mal[columns_to_keep]
    else:
        mal = mal[
            [
                "CASES",
                "DISTRICT",
                "AGE",
                "LON",
                "LAT",
                "ALT",
                "POPAGE",
                "TOTALPOP",
                "CALMONTH",
                "CALYEAR",
                "PRECIP",
                "TEMP",
            ]
        ]

    mal["AGE"] = "ALL"

    # get year range
    mal = mal[(mal["CALYEAR"] >= year_start) & (mal["CALYEAR"] <= year_end)]

    # remove district 10 for lack of data/trust
    mal = mal[~mal["DISTRICT"].isin(dist_to_remove)]

    # get malawi case and demographic info
    # initialize nan arrays and then populate cases/pop with the data
    num_years = year_end - year_start + 1
    yr = np.arange(year_start, year_end + 1, 1 / 12)

    # create datetime object based on calmonth/calyear columns. note: using day=15 for all months as a dummy day
    mal["DATETIME"] = pd.to_datetime(
        mal["CALYEAR"].astype(str) + "-" + mal["CALMONTH"].astype(str) + "-15"
    )

    # get cases rate per 1000 pop
    mal["CASERATE1000"] = (
        1e3 * mal["CASES"] / mal["TOTALPOP"]
    )  # case rate per 1000 people

    return mal
