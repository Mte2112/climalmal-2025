import pandas as pd
import geopandas as gpd
from scipy.signal import butter, lfilter, freqz, filtfilt
import numpy as np
import xarray as xr
import rasterio
from shapely.geometry import mapping
from sklearn.linear_model import LinearRegression

# modify filter parameters based on period input
def calculate_cutoff(per, fs):
    return 1 / per * fs


# define the filter function
def butter_filter(data, cutoff, fs, order=5, btype="low"):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype=btype)
    return filtfilt(b, a, data)


def pad_df(data):
    ns = len(data) // 2
    sim1 = data.iloc[:ns].iloc[::-1]  # flipped version of the first half
    sim2 = data.iloc[-ns:].iloc[::-1]  # flipped version of the last half
    return pd.concat([sim1, data, sim2]).reset_index(drop=True)


def pad_array(arr):
    # determine the number of samples to pad (half the length of the array)
    ns = len(arr) // 2

    # create reversed versions of the first and last halves
    sim1 = arr[:ns][::-1]  # reverse the first half
    sim2 = arr[-ns:][::-1]  # reverse the last half

    # concatenate the padded data
    padded_array = np.concatenate([sim1, arr, sim2])
    return padded_array


def pad_xr(data):
    # get the number of samples to pad from the time dimension
    ns = data.sizes["time"] // 2  # number of samples to pad
    # create the flipped versions for padding
    sim1 = data.isel(time=slice(0, ns)).isel(
        time=slice(None, None, -1)
    )  # flipped version of the first half
    sim2 = data.isel(time=slice(-ns, None)).isel(
        time=slice(None, None, -1)
    )  # flipped version of the last half
    # concatenate along the time dimension
    return xr.concat([sim1, data, sim2], dim="time")


def apply_filter_cases(group, casevar, filter_type="low", per=24, fs=1, order=3):
    cutoff = 1 / per * fs  # calculate cutoff
    padded_data = pad_df(group[casevar])
    smoothed_data = butter_filter(padded_data.values, cutoff, fs, order, filter_type)
    # group['CASEANOM_smoothed'] = smoothed_data[len(group) // 2:len(group) // 2 + len(group)]
    ns = len(group) // 2  # number of samples padded
    group["CASES_smoothed"] = smoothed_data[
        ns : ns + len(group)
    ]  # remove the padded data
    return group


def area_weighted_mean(data_array, lat_bounds, lon_bounds):
    """
    calculate the area-weighted mean for a given SST data subset.

    Parameters:
    sst_data (xarray.DataArray): the SST data to analyze.
    lat_bounds (tuple): latitude bounds for selection (min_lat, max_lat).
    lon_bounds (tuple): longitude bounds for selection (min_lon, max_lon).

    Returns:
    xarray.DataArray: the area-weighted mean.
    """
    if "latitude" in data_array.dims:
        data_array = data_array.rename({"latitude": "lat"})
        data_array.lat.attrs["standard_name"] = "lat"
        data_array.lat.attrs["long_name"] = "lat"
    if "longitude" in data_array.dims:
        data_array = data_array.rename({"longitude": "lon"})
        data_array.lon.attrs["standard_name"] = "lon"
        data_array.lon.attrs["long_name"] = "lon"

    # select the region based on latitude and longitude bounds
    data_region = data_array.where(
        (data_array.lat >= lat_bounds[0])
        & (data_array.lat <= lat_bounds[1])
        & (data_array.lon >= lon_bounds[0])
        & (data_array.lon <= lon_bounds[1]),
        drop=True,
    )

    # compute the cosine of latitudes for area weighting
    weights = np.cos(np.deg2rad(data_region["lat"]))
    # normalize the weights so they sum to 1
    weights /= weights.sum()
    # apply the weights along the latitude dimension and compute the weighted mean
    weighted_mean = data_region.weighted(weights).mean(("lon", "lat"))

    return weighted_mean


def smooth_da_ts(da, namevar, per=(5 * 12), fs=1, order=3):
    """
    low-pass filter data array timeseries

    Returns:
    xarray.Dataset: filtered ds
    """

    cutoff = 1 / per * fs
    da_val = da.values
    # pad data for smoothing
    da_padded = pad_array(da_val)
    # apply the Butterworth filter to the dmi data
    da_smoothed = butter_filter(da_padded, cutoff, fs, order, btype="low")
    # now drop the padded data
    ns = len(da_val) // 2  # number of samples padded
    da_smoothed = da_smoothed[ns : ns + len(da_val)]
    smoothed_ds = xr.Dataset(
        data_vars={namevar: (["time"], da_smoothed)}, coords={"time": da.time}
    )

    return smoothed_ds


def clip(ds, shape, buffer=0.1, londim="lon", latdim="lat", drop=True):
    """
    clip an xarray dataset within the bounds of a shape
    """

    # use the shapefile to clip the dataset, buffering the boundaries
    shape_buffered = gpd.GeoDataFrame(geometry=shape.buffer(buffer)).to_crs("EPSG:4326")

    ds.rio.set_spatial_dims(x_dim=londim, y_dim=latdim, inplace=True)
    ds.rio.write_crs("EPSG:4326", inplace=True)
    ds_clipped = ds.rio.clip(
        shape_buffered.geometry.apply(mapping), shape_buffered.crs, drop=drop
    )

    return ds_clipped


def detrend_dim(da, dim, deg=1):
    """detrend data"""
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(da[dim], p.polyfit_coefficients)
    return da - fit


def detrend(da, dims, deg=1):
    """detrend along multiple dimensions"""
    # only valid for linear detrending (deg=1)
    da_detrended = da
    for dim in dims:
        da_detrended = detrend_dim(da_detrended, dim, deg=deg)
    return da_detrended


def normalize(ts):
    """normalize timeseries to unit standard deviation"""
    return (ts - np.mean(ts)) / np.std(ts)


def regmap_ts_simple_fast(predictand, predictor):
    """dask-leveraged regression of the spatiotemporal field onto a timeseries"""
    # check if numpy array
    if isinstance(predictor, np.ndarray):
        ts = predictor
    else:
        ts = predictor.values

    # normalize time series to get in terms of standard deviation rather than nominal value
    ts = normalize(ts)

    # add intercept/constant
    ts = np.stack((np.ones_like(ts), ts), axis=-1)  # shape: [time, 2]

    # switch lon/lat to long name for consistency below
    if "lat" in predictor.dims:
        predictor = predictor.rename({"lat": "latitude"})
        predictor = predictor.rename({"lon": "longitude"})
    if "lat" in predictand.dims:
        predictand = predictand.rename({"lat": "latitude"})
        predictand = predictand.rename({"lon": "longitude"})
    # add predictor dimension
    predictor_broadcasted = predictor.expand_dims(
        {"latitude": predictand.latitude, "longitude": predictand.longitude}
    ).transpose("time", "latitude", "longitude")

    # chunk predictand for dask
    pdctd_chunk = predictand.chunk({"latitude": -1, "longitude": -1})

    # apply regression at each grid point func
    # note to self: quick and dirty for now. move out later
    def apply_regression(grid_point_ts, ts):
        ts = ts.reshape(-1, 1)  # check if 2D

        # does ts have more than one feature?
        if ts.shape[1] == 1:
            # if only one feature, the model will return only one coefficient
            model = LinearRegression().fit(ts, grid_point_ts)
            return model.coef_[0]
        else:
            # if multiple features, handle it as usual
            model = LinearRegression().fit(ts, grid_point_ts)
            return model.coef_[1]

    # apply regression
    regmap = xr.apply_ufunc(
        apply_regression,
        pdctd_chunk,
        predictor_broadcasted,
        input_core_dims=[["time"], ["time"]],
        output_core_dims=[[]],  # output is scalar (slope)
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64],
    )

    return regmap

def apply_filter_xr(da, per, fs, cutoff, order):
    '''apply LPF to xarray data array'''
    
    # apply low-pass filter to SST data
    padded_data = pad_xr(da) # pad data to manage tails 
    
    # apply the filter to each grid cell. this may take a couple minutes
    filtered_data = xr.apply_ufunc(
        butter_filter,
        padded_data,
        input_core_dims=[['time']],  # filter along time
        output_core_dims=[['time']],
        vectorize=True,
        dask='allowed', 
        kwargs={'cutoff': cutoff, 'fs': fs, 'order': order}
    )
    # drop padded data 
    ns = da.sizes['time'] // 2  # number of samples padded
    filtered_data = filtered_data.isel(time=slice(ns, ns + da.sizes['time'])) 

    return filtered_data


# get district correlations
def get_dist_corr(clim_da, var, dist_shape, mal):
    '''get correlation of district malaria climate variable within that district'''
    
    corrs = []
    for dnum in dist_shape["DISTRICT"].unique():
        district = dist_shape[dist_shape["DISTRICT"] == dnum]
        # use the malawi shapefile to mask the dataset, buffering the boundaries by about 10 km to capture overlapping grid cells on boundaries
        dist_buffered = gpd.GeoDataFrame(geometry=district.buffer(0.1)).to_crs("EPSG:4326")
    
        # clip climate var within district boundaries
        if "lon" in clim_da.dims:
            lon = "lon"
            lat = "lat"
        else:
            lon = "longitude"
            lat = "latitude"
        clim_da.rio.set_spatial_dims(x_dim=lon, y_dim=lat, inplace=True)
        clim_da.rio.write_crs("EPSG:4326", inplace=True)
        clim_da_dist = clim_da.rio.clip(dist_buffered.geometry.apply(mapping), dist_buffered.crs, drop=True)
    
        # smooth clim_da 
        clim_da_dist_mean = area_weighted_mean(clim_da_dist[var], (-90, 90), (0, 360))
        clim_da_dist_mean_smooth = smooth_da_ts(clim_da_dist_mean, var, per=(3*12))
    
        # smooth cases in district
        dist_cases = mal[mal["DISTRICT"] == dnum].groupby("DATETIME")["CASES"].sum()
        da_cases_dist = xr.DataArray(dist_cases.values)
        cases_dist_ds = xr.Dataset(
                data_vars={"cases": (["time"], da_cases_dist.values)},
                coords={"time": dist_cases.index.values}
            )
        
        da_cases_dist_smooth = smooth_da_ts(cases_dist_ds.cases, "cases", per=(2*12))
        da_cases_dist_smooth["time"] = clim_da_dist_mean_smooth.time
    
        # get correlation between district cases and average clim_da 
        corr = xr.corr(da_cases_dist_smooth.cases, clim_da_dist_mean_smooth[var])
        corrs.append(corr.values.item())

    return corrs