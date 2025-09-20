"""
Computations to produce the core results in Elling et al., 2025.
Absent are the CMIP6 analyses which are included in a separate script
Public climate datasets will need to be downloaded via the links below and paths set directly below
Note: there are a lot of calculations and a lot of data, so it may take the script a while to run
but this script is intended to be run only once, which is the purpose of writing (hopefully not a painful amount of) intermediate outputs
Preprocessed data can then be imported from preprocessed/ w/o rerunning
See README for information
"""

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# NOTE: you must modify vars below to match your paths #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# malaria data - see Zenodo
MALARIA_DATA_PATH = ""

# shapefiles - provided in data/shapefiles. Unzip if compressed
MALAWI_SHAPE_PATH = "data/shapefiles/mwi_admbnda_adm0_nso_hotosm_20230405.shp"
DISTRICT_SHAPE_PATH = "data/shapefiles/maldists/MWI_adm1.shp"
LAKE_MALAWI_PATH = "data/shapefiles/lake_malawi.shp"

# topography - publically available https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/tiled/
ETOPO_PATH = ""

# climate data paths - see link 4 public download
SST_PATH = ""  # https://climatedataguide.ucar.edu/climate-data/sst-data-noaa-optimal-interpolation-oi-sst-analysis-version-2-oisstv2-1x1
# era5 data https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview
WIND_PATH = ""
GEOP_300_PATH = ""
GEOP_850_PATH = ""
VIMFD_PATH = (
    ""
)
TEMP_PATH = ""
PRECIP_PATH = "."  # https://www.chc.ucsb.edu/data/chirps
GRACE_PATH = ""  # https://disc.gsfc.nasa.gov/datasets/GRACEDADM_CLSM025GL_7D_3.0/summary?keywords=grace%20soil%20moisture
GRACE_MALAWI_PATH = ""  # preprocessed from dataset above. Malawi clipped. see analysis_tools.py function clip()
SOIL_TEXTURE_PATH = ""  # supplementary https://ldas.gsfc.nasa.gov/gldas/soils

# output preprocessed data to this directory
OUTPUT_DIR = "preprocessed/"

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#                  End path setting                    #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

import analysis_tools as atools
import format_data as fdata
import geopandas as gpd
import matplotlib as mpl
import numpy as np
import os
import pandas as pd
import warnings
import xarray as xrs
from eofs.standard import Eof
from eofs.xarray import Eof as xEof
from matplotlib.colors import BoundaryNorm, ListedColormap
from scipy.stats import zscore
from shapely.geometry import mapping

# uncomment if parallelization wanted, but you will prob need to chunk the xr datasets and do some restructuring
# from concurrent.futures import ProcessPoolExecutor
# import dask
# dask.config.set(scheduler='threads', num_workers=4)

# mpl.rcParams['font.family'] = 'serif'
warnings.filterwarnings("ignore")


def create_output_directory(output_dir):
    """create output directory for preprocessed data"""
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def load_malaria_data(malaria_path):
    """load and process malaria incidence data"""
    print("Processing malaria data...")

    # date range for analysis
    y1 = 2005
    y2 = 2023
    casevar = "CASES"

    # read and format malaria data
    # note: district 10 is dropped for insufficient data
    mal = fdata.process_malaria_data(
        malaria_path,
        y1,
        y2,
        2,
        [10],
        ["CASES", "DISTRICT", "LON", "LAT", "TOTALPOP", "CALMONTH", "CALYEAR"],
    )

    return mal, y1, y2, casevar


def load_supporting_data(malawi_shape_path, district_shape_path, lake_path, etopo_path):
    """load shape files and topographic data for map construction"""
    print("Loading supporting geographic data...")

    # read in Malawi shape file, whole country
    malawi_shape = gpd.read_file(malawi_shape_path)

    # read Malawi district boundaries
    dist_shape = gpd.read_file(district_shape_path)
    dist_shape = dist_shape[dist_shape["ID_1"] != 10].rename(
        columns={"ID_1": "DISTRICT"}
    )

    # lake malawi shape file - for maps
    lake_malawi = gpd.read_file(lake_path)

    # read in topo (elevation)
    etopo01 = xr.open_dataset(etopo_path)

    # subset the data around Malawi to speed up plotting performance
    lon_mask = (etopo01.x >= 10) & (etopo01.x <= 60)
    lat_mask = (etopo01.y >= -30) & (etopo01.y <= 0)
    etopo01 = etopo01.where((lon_mask) & (lat_mask))

    return malawi_shape, dist_shape, lake_malawi, etopo01


def compute_malaria_statistics(mal, dist_shape, casevar):
    """compute national and district-level malaria statistics including EOFs"""
    print("Computing malaria statistics and EOFs...")

    # get case rate per 1000 people
    malavg = mal.groupby("DISTRICT").mean("CASERATE1000")

    # add district polygon to malaria df
    malavg = malavg.reset_index()
    malavg["geometry"] = pd.merge(malavg, dist_shape, how="left", on="DISTRICT")[
        "geometry"
    ]
    malavg = gpd.GeoDataFrame(malavg)
    malavgscale = malavg.copy()
    malavgscale["CASES"] = malavgscale["CASES"] * 1e-4

    # get district level cases and apply a low-pass filter
    mal = (
        mal.groupby("DISTRICT")
        .apply(
            lambda g: atools.apply_filter_cases(g, casevar, filter_type="low", per=24)
        )
        .reset_index(drop=True)
    )

    # national total monthly cases
    nattot = mal.groupby("DATETIME")["CASES"].sum()

    # get malaria EOFs
    cases_pivot = mal.pivot_table(index="DATETIME", columns="DISTRICT", values="CASES")
    cases_array = cases_pivot.values

    # do eof analysis, getting first 5 eofs
    eof_solver = Eof(cases_array)
    eofs = eof_solver.eofs(neofs=5)
    explained_variance = eof_solver.varianceFraction(neigs=5) * 100

    # add eofs to df
    malavg["EOF1"] = eofs[0]
    malavg["EOF2"] = eofs[1]
    malavg["EOF3"] = eofs[2]
    malavg["EOF4"] = eofs[3]
    malavg["EOF5"] = eofs[4]

    # get principal components (5)
    pcs = eof_solver.pcs(npcs=5)
    pc1 = pcs[:, 0]
    pc2 = pcs[:, 1]

    # low-pass filter malaria PC's
    pc1_sm = atools.apply_filter_cases(
        pd.DataFrame({"CASES": pc1}), "CASES", filter_type="low", per=24
    )
    pc2_sm = atools.apply_filter_cases(
        pd.DataFrame({"CASES": pc2}), "CASES", filter_type="low", per=24
    )

    # set pc1 and pc2 variable values to lpf data because lpf pcs exclusively used in figures
    pc1 = pc1_sm["CASES_smoothed"].values
    pc2 = pc2_sm["CASES_smoothed"].values

    # get EOFs with 24mo smooth
    pcs_cases_24 = eof_solver.pcs(npcs=5)
    pc1_cases_24 = pcs_cases_24[:, 0]
    pc2_cases_24 = pcs_cases_24[:, 1]
    sm_pc1_cases_24 = (
        pd.DataFrame({"DISTRICT": 1, casevar: pc1_cases_24})
        .groupby("DISTRICT")
        .apply(
            lambda g: atools.apply_filter_cases(g, casevar, filter_type="low", per=24)
        )
        .reset_index(drop=True)["CASES_smoothed"]
    )
    sm_pc2_cases_24 = (
        pd.DataFrame({"DISTRICT": 1, casevar: pc2_cases_24})
        .groupby("DISTRICT")
        .apply(
            lambda g: atools.apply_filter_cases(g, casevar, filter_type="low", per=24)
        )
        .reset_index(drop=True)["CASES_smoothed"]
    )

    print(f"EOFs explained variance 1-5: {explained_variance}")
    print(f"sum explained variance eof1+eof2: {sum(explained_variance[0:2])}")

    return (
        mal,
        malavg,
        malavgscale,
        nattot,
        cases_pivot,
        pc1,
        pc2,
        sm_pc1_cases_24,
        sm_pc2_cases_24,
        explained_variance,
    )


def process_sst_data(sst_path, y1, y2, pc1, pc2):
    """process sea surface temperature data and compute correlations"""
    print("Processing SST data...")

    # read in NOAA OI SST dataset
    ds_sst = xr.open_dataset(sst_path).sel(time=slice(str(y1), str(y2)))
    ds_sst = fdata.fliplon(ds_sst)

    # we previously identified that TA and IO oceans are our main drivers
    # save these timeseries before low-pass filtering
    tropatl_weighted_nosmooth = atools.area_weighted_mean(ds_sst, (-8, 2), (-40, -10))
    io_sst_weighted_nosmooth = atools.area_weighted_mean(ds_sst, (-26, 2), (64, 92))

    # get unfiltered SST anomalies
    io_base_nosmooth = io_sst_weighted_nosmooth.groupby("time.month").mean()
    tropatl_base_nosmooth = tropatl_weighted_nosmooth.groupby("time.month").mean()
    io_anom_nosmooth = io_sst_weighted_nosmooth.groupby("time.month") - io_base_nosmooth
    tropatl_anom_nosmooth = (
        tropatl_weighted_nosmooth.groupby("time.month") - tropatl_base_nosmooth
    )

    # apply low-pass filter to SST data
    per = 5 * 12
    fs = 1
    cutoff = 1 / per
    order = 3
    filtered_sst = atools.apply_filter_xr(ds_sst.sst, per, fs, cutoff, order)

    # get the mean of the lpf TA and IO regions
    tropatl_weighted_smooth = atools.area_weighted_mean(
        filtered_sst, (-8, 2), (-40, -10)
    )
    io_weighted_smooth = atools.area_weighted_mean(filtered_sst, (-26, 2), (64, 92))

    # get correlation sst and pc1 malaria
    pc1_da = xr.DataArray(pc1, dims="time", coords={"time": filtered_sst.time})
    correlation_map_pc1 = xr.corr(filtered_sst, pc1_da, dim="time")

    # same for pc2
    pc2_da = xr.DataArray(pc2, dims="time", coords={"time": filtered_sst.time})
    correlation_map_pc2 = xr.corr(filtered_sst, pc2_da, dim="time")

    return (
        tropatl_weighted_nosmooth,
        io_sst_weighted_nosmooth,
        tropatl_anom_nosmooth,
        io_anom_nosmooth,
        tropatl_weighted_smooth,
        io_weighted_smooth,
        correlation_map_pc1,
        correlation_map_pc2,
    )


def process_wind_data(wind_path, y1, y2, tropatl_anom_nosmooth, io_anom_nosmooth):
    """process ERA5 wind data and compute regressions"""
    print("Processing wind data...")

    # read in era5 wind data and format
    wind850mb = (
        xr.open_mfdataset(wind_path).sortby("time").sel(time=slice(str(y1), str(y2)))
    )
    wind850mb = fdata.fliplon(wind850mb)
    wind850mb = wind850mb.sortby(["time", "longitude", "latitude"])

    # u-v vector formatting
    u = wind850mb["u"]
    v = wind850mb["v"]
    u_mean = u.mean("time")
    v_mean = v.mean("time")
    u_anom = u.groupby("time.month") - u.groupby("time.month").mean()
    v_anom = v.groupby("time.month") - v.groupby("time.month").mean()
    u_anom = u_anom.chunk({"latitude": -1, "longitude": -1, "time": -1})
    v_anom = v_anom.chunk({"latitude": -1, "longitude": -1, "time": -1})

    # run wind regression, SST as predictor
    u_reg_tropatl = atools.regmap_ts_simple_fast(
        u_anom, tropatl_anom_nosmooth.sst
    ).compute()
    v_reg_tropatl = atools.regmap_ts_simple_fast(
        v_anom, tropatl_anom_nosmooth.sst
    ).compute()
    u_reg_io = atools.regmap_ts_simple_fast(u_anom, io_anom_nosmooth.sst).compute()
    v_reg_io = atools.regmap_ts_simple_fast(v_anom, io_anom_nosmooth.sst).compute()

    return u, v, u_mean, v_mean, u_reg_tropatl, v_reg_tropatl, u_reg_io, v_reg_io


def process_geopotential_data(
    geop_300_path,
    geop_850_path,
    y1,
    y2,
    tropatl_anom_nosmooth,
    io_anom_nosmooth,
    tropatl_weighted_smooth,
    tropatl_weighted_nosmooth,
):
    """process geopotential height data"""
    print("Processing geopotential height data...")

    # 300 mb geopotential height
    geop_300 = xr.open_dataset(geop_300_path).sel(time=slice(str(y1), str(y2)))
    geop_300 = geop_300 / 9.81  # get geopotential height from geopotential
    geop_300["time"] = tropatl_weighted_smooth.time  # align time coords
    geop_300_anom = (
        geop_300.groupby("time.month") - geop_300.groupby("time.month").mean()
    )

    # 850 mb geopotential
    geop_850 = xr.open_dataset(geop_850_path).sel(time=slice(str(y1), str(y2)))
    geop_850 = geop_850 / 9.81  # get geopotential height
    geop_850["time"] = tropatl_weighted_nosmooth.time
    geop_850_height_mean = geop_850.z.mean("time")
    geop_850_anom = (
        geop_850.groupby("time.month") - geop_850.groupby("time.month").mean()
    )

    # regress geop on sst
    geop_300_anom = geop_300_anom.chunk({"time": -1})
    reg_geop_300_anom_tropatl = atools.regmap_ts_simple_fast(
        geop_300_anom.z, tropatl_anom_nosmooth.sst
    ).compute()
    reg_geop_300_anom_io = atools.regmap_ts_simple_fast(
        geop_300_anom.z, io_anom_nosmooth.sst
    ).compute()

    geop_850_anom = geop_850_anom.chunk({"time": -1})
    reg_geop_850_anom_tropatl = atools.regmap_ts_simple_fast(
        geop_850_anom.z, tropatl_anom_nosmooth.sst
    ).compute()
    reg_geop_850_anom_io = atools.regmap_ts_simple_fast(
        geop_850_anom.z, io_anom_nosmooth.sst
    ).compute()

    return (
        geop_850_height_mean,
        reg_geop_300_anom_tropatl,
        reg_geop_300_anom_io,
        reg_geop_850_anom_tropatl,
        reg_geop_850_anom_io,
    )


def process_precipitation_data(
    precip_path,
    y1,
    y2,
    malawi_shape,
    tropatl_weighted_smooth,
    io_weighted_smooth,
    dist_shape,
    mal,
):
    """process CHIRPS precipitation data"""
    print("Processing precipitation data...")

    # read chirps precip
    precip_raw = (
        xr.open_mfdataset(precip_path)
        .sortby("time")
        .sel(
            longitude=slice(25, 45),
            latitude=slice(-20, -8),
            time=slice(str(y1), str(y2)),
        )
        .compute()
    )
    precip = precip_raw.copy()

    # get climatology stats
    malawi_shape_buffered = gpd.GeoDataFrame(geometry=malawi_shape.buffer(0.01)).to_crs(
        "EPSG:4326"
    )
    precip_mal = precip.copy()
    precip_mal.groupby("time.month")
    precip_mal.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
    precip_mal.rio.write_crs("EPSG:4326", inplace=True)
    precip_mal = precip_mal.rio.clip(
        malawi_shape_buffered.geometry.apply(mapping),
        malawi_shape_buffered.crs,
        drop=True,
    )

    # get average precip
    precip_mal = atools.area_weighted_mean(precip_mal["precip"], (-90, 90), (0, 360))
    precip_mal_nov_mar = precip_mal.where(
        precip_mal.time.dt.month.isin([11, 12, 1, 2, 3]), drop=True
    )
    precip_mal_jun_sep = precip_mal.where(
        precip_mal.time.dt.month.isin([6, 7, 8, 9]), drop=True
    )

    print(f"wet season precip: {precip_mal_nov_mar.mean('time').values}")
    print(f"dry season precip: {precip_mal_jun_sep.mean('time').values}")

    # format precip data for dask-leveraged regression
    precip = precip["precip"]
    precip = precip.where(precip != precip.encoding["missing_value"], np.nan)
    precip = precip.fillna(-99999)
    mask = precip != -99999
    nan_mask = mask.all(dim="time")
    tropatl_weighted_smooth["time"] = precip["time"]
    tropatl_weighted_smooth = tropatl_weighted_smooth.chunk({"time": -1})
    precip = precip.chunk({"time": -1})

    # run regression
    precip_reg_TA = atools.regmap_ts_simple_fast(
        precip, tropatl_weighted_smooth
    ).compute()
    precip_reg_IO = atools.regmap_ts_simple_fast(precip, io_weighted_smooth).compute()

    # mask all data gaps
    precip_reg_masked_TA = precip_reg_TA.where(nan_mask)
    precip_reg_masked_IO = precip_reg_IO.where(nan_mask)

    # precip dist corr
    precip = precip_raw.copy()
    corrs = atools.get_dist_corr(precip, "precip", dist_shape, mal)
    dist_shape["corr_precip"] = corrs
    malavg_precip = (
        mal.reset_index()
    )  # this should use malavg but keeping consistent with original
    malavg_precip["corr_precip"] = pd.merge(
        malavg_precip, dist_shape, how="left", on="DISTRICT"
    )["corr_precip"]

    return precip_reg_masked_TA, precip_reg_masked_IO, malavg_precip


def process_vimfd_data(vimfd_path, y1, y2, tropatl_anom_nosmooth, io_anom_nosmooth):
    """process vertically integrated moisture flux divergence"""
    print("Processing VIMFD data...")

    # read and format vimfd
    vimdf = xr.open_dataset(vimfd_path).sel(time=slice(str(y1), str(y2)))

    # get anomaly
    vimdf_anom = vimdf.groupby("time.month") - vimdf.groupby("time.month").mean()

    # set up and run vimfd on TA/IO regression
    vimdf_anom = vimdf_anom.chunk({"time": -1})
    reg_vimdf_TA = atools.regmap_ts_simple_fast(
        vimdf_anom.vimdf, tropatl_anom_nosmooth.sst
    ).compute()
    reg_vimdf_IO = atools.regmap_ts_simple_fast(
        vimdf_anom.vimdf, io_anom_nosmooth.sst
    ).compute()

    # spatially smooth data for plotting
    new_lon = np.linspace(
        reg_vimdf_TA.longitude[0].item(),
        reg_vimdf_TA.longitude[-1].item(),
        reg_vimdf_TA.sizes["longitude"] * 4,
    )
    new_lat = np.linspace(
        reg_vimdf_TA.latitude[0].item(),
        reg_vimdf_TA.latitude[-1].item(),
        reg_vimdf_TA.sizes["latitude"] * 4,
    )
    reg_vimdf_TA_smoothed = reg_vimdf_TA.interp(latitude=new_lat, longitude=new_lon)
    reg_vimdf_IO_smoothed = reg_vimdf_IO.interp(latitude=new_lat, longitude=new_lon)

    return reg_vimdf_TA_smoothed, reg_vimdf_IO_smoothed


def process_soil_moisture_data(
    grace_path,
    grace_malawi_path,
    y1,
    y2,
    tropatl_anom_nosmooth,
    io_anom_nosmooth,
    sm_pc1_cases_24,
    sm_pc2_cases_24,
    casevar,
    dist_shape,
    mal,
):
    """process GRACE soil moisture data"""
    print("Processing soil moisture data...")

    # read GRACE surface soil moisture
    gracevar = "sfsm_inst"
    soil_moisture_malreg = fdata.process_nc_malawi(
        grace_path, y1, y2, [25, 45, -20, -8], "mean"
    ).sortby("time")
    soilm_mean_region = atools.area_weighted_mean(
        soil_moisture_malreg[gracevar], (-90, 90), (0, 360)
    )

    # malawi only grace data (preprocessed for speed)
    soil_moisture = xr.open_dataset(grace_malawi_path).sel(time=slice(str(y1), str(y2)))
    soilm_mean = atools.area_weighted_mean(soil_moisture[gracevar], (-90, 90), (0, 360))

    # get anomalies & regression on SST anomalies
    soilm_anom = (
        soil_moisture_malreg.groupby("time.month")
        - soilm_mean_region.groupby("time.month").mean()
    )
    tropatl_weighted_nosmooth_anom = tropatl_anom_nosmooth  # already computed
    io_sst_weighted_nosmooth_anom = io_anom_nosmooth  # already computed

    # set up regression
    soilm_anom = soilm_anom.chunk({"lat": -1, "lon": -1, "time": -1})[
        "sfsm_inst"
    ].rename({"lat": "latitude", "lon": "longitude"})
    # handle nans
    soilm_anom = soilm_anom.fillna(-99999)
    mask = soilm_anom != -99999
    nan_mask = mask.all(dim="time")

    tropatl_weighted_nosmooth_anom["time"] = soilm_anom["time"]
    tropatl_weighted_nosmooth_anom = tropatl_weighted_nosmooth_anom.chunk({"time": -1})
    io_sst_weighted_nosmooth_anom["time"] = soilm_anom["time"]
    io_sst_weighted_nosmooth_anom = io_sst_weighted_nosmooth_anom.chunk({"time": -1})

    # soil regressionsz
    soil_reg_tropatl = atools.regmap_ts_simple_fast(
        soilm_anom, tropatl_weighted_nosmooth_anom.sst
    ).compute()
    soil_reg_io = atools.regmap_ts_simple_fast(
        soilm_anom, io_sst_weighted_nosmooth_anom.sst
    ).compute()
    soil_reg_tropatl_masked = soil_reg_tropatl.where(nan_mask)
    soil_reg_io_masked = soil_reg_io.where(nan_mask)

    # now set up soil moisture-malaria comparisons
    # get 36 mo smoothed soil moisture
    padded_grace = atools.pad_xr(soil_moisture[gracevar])
    per = 3 * 12
    fs = 1
    cutoff = 1 / per
    order = 3
    filtered_grace = atools.apply_filter_xr(
        soil_moisture[gracevar], per, fs, cutoff, order
    )

    # get Malawi average
    soilm_smth36mo_avg = atools.area_weighted_mean(filtered_grace, (-90, 90), (0, 360))

    # get correlation soil moisture 36 mo smoothed and pc1 24 mo smoothed
    sm_pc1_cases_24_da = xr.DataArray(
        sm_pc1_cases_24, dims="time", coords={"time": filtered_grace.time}
    )
    correlation_map_soilm_pc1 = xr.corr(
        filtered_grace, sm_pc1_cases_24_da, dim="time"
    ).compute()

    # get correlation soil and pc2 malaria
    sm_pc2_cases_24_da = xr.DataArray(
        sm_pc2_cases_24, dims="time", coords={"time": filtered_grace.time}
    )
    correlation_map_soilm_pc2 = xr.corr(
        filtered_grace, sm_pc2_cases_24_da, dim="time"
    ).compute()

    # soil-malaria correlation at district level
    corrs = atools.get_dist_corr(soil_moisture, "sfsm_inst", dist_shape, mal)
    dist_shape["corr_soilm"] = corrs
    malavg_soilm = mal.reset_index()  # keeping consistent with original structure
    malavg_soilm["corr_soilm"] = pd.merge(
        malavg_soilm, dist_shape, how="left", on="DISTRICT"
    )["corr_soilm"]

    # get soilm pcs for PC comparisons
    soilm_detrended = soil_moisture.sfsm_inst - soil_moisture.sfsm_inst.mean("time")
    soilm_detrended = atools.detrend(soilm_detrended, ["time"], deg=1)
    coslat = np.cos(np.deg2rad(soilm_detrended.coords["lat"].values))
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = xEof(soilm_detrended, weights=wgts)
    pc1_soilm_mal = solver.pcs(npcs=5).sel(mode=0).values
    pc2_soilm_mal = solver.pcs(npcs=5).sel(mode=1).values
    print(f"soil moisture eofs EV: {solver.varianceFraction(neigs=5) * 100}")

    # smoothing after eof
    sm_soilm_pc1 = (
        pd.DataFrame({"DISTRICT": 1, casevar: pc1_soilm_mal})
        .groupby("DISTRICT")
        .apply(
            lambda g: atools.apply_filter_cases(g, casevar, filter_type="low", per=36)
        )
        .reset_index(drop=True)["CASES_smoothed"]
    )
    sm_soilm_pc2 = (
        pd.DataFrame({"DISTRICT": 1, casevar: pc2_soilm_mal})
        .groupby("DISTRICT")
        .apply(
            lambda g: atools.apply_filter_cases(g, casevar, filter_type="low", per=36)
        )
        .reset_index(drop=True)["CASES_smoothed"]
    )

    return (
        soil_reg_tropatl_masked,
        soil_reg_io_masked,
        correlation_map_soilm_pc1,
        correlation_map_soilm_pc2,
        malavg_soilm,
        sm_soilm_pc1,
        sm_soilm_pc2,
    )


def process_temperature_data(
    temp_path,
    y1,
    y2,
    malawi_shape,
    tropatl_anom_nosmooth,
    io_anom_nosmooth,
    dist_shape,
    mal,
    io_weighted_smooth,
):
    """process ERA5 2m temperature data"""
    print("Processing temperature data...")

    # read in dataset
    temp_global = xr.open_dataset(temp_path).sel(time=slice(str(y1), str(y2)))
    temp = temp_global.where(
        (temp_global.longitude >= 25)
        & (temp_global.longitude <= 45)
        & (temp_global.latitude >= -20)
        & (temp_global.latitude <= -8)
    ).drop(["number", "step", "surface", "valid_time"])

    # get climatology stats
    malawi_shape_buffered = gpd.GeoDataFrame(geometry=malawi_shape.buffer(0.01)).to_crs(
        "EPSG:4326"
    )
    temp_mal = temp.copy()
    temp_mal.groupby("time.month")
    temp_mal.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
    temp_mal.rio.write_crs("EPSG:4326", inplace=True)
    temp_mal = temp_mal.rio.clip(
        malawi_shape_buffered.geometry.apply(mapping),
        malawi_shape_buffered.crs,
        drop=True,
    )

    # get average temp
    temp_mal = atools.area_weighted_mean(temp_mal["t2m"], (-90, 90), (0, 360))
    temp_mal_nov = temp_mal.where(temp_mal.time.dt.month == 11, drop=True)
    temp_mal_jul = temp_mal.where(temp_mal.time.dt.month == 7, drop=True)
    print(f"nov temp: {temp_mal_nov.mean('time').values - 273.15}")
    print(f"july temp: {temp_mal_jul.mean('time').values - 273.15}")

    # set up and run temp-sst regression
    temp_anom = temp.groupby("time.month") - temp.groupby("time.month").mean()
    temp_anom = temp_anom.chunk({"latitude": -1, "longitude": -1, "time": -1})["t2m"]

    # handle nans
    temp_anom = temp_anom.fillna(-99999)
    mask = temp_anom != -99999
    nan_mask = mask.all(dim="time")

    tropatl_weighted_nosmooth_anom = tropatl_anom_nosmooth
    tropatl_weighted_nosmooth_anom["time"] = temp_anom["time"]
    tropatl_weighted_nosmooth_anom = tropatl_weighted_nosmooth_anom.chunk({"time": -1})
    io_sst_weighted_nosmooth_anom = io_anom_nosmooth
    io_sst_weighted_nosmooth_anom["time"] = temp_anom["time"]
    io_sst_weighted_nosmooth_anom = io_sst_weighted_nosmooth_anom.chunk({"time": -1})

    # temp regressions
    temp_reg_tropatl = atools.regmap_ts_simple_fast(
        temp_anom, tropatl_weighted_nosmooth_anom.sst
    ).compute()
    temp_reg_io = atools.regmap_ts_simple_fast(
        temp_anom, io_sst_weighted_nosmooth_anom.sst
    ).compute()
    temp_reg_tropatl_masked = temp_reg_tropatl.where(nan_mask)
    temp_reg_io_masked = temp_reg_io.where(nan_mask)

    # temp dist corr
    corrs = atools.get_dist_corr(temp, "t2m", dist_shape, mal)
    dist_shape["corr_temp"] = corrs
    malavg_temp = mal.reset_index()  # keeping consistent with original structure
    malavg_temp["corr_temp"] = pd.merge(
        malavg_temp, dist_shape, how="left", on="DISTRICT"
    )["corr_temp"]

    # supplementary figure - temperature advection
    # shows that IO is warmer than malawi when positive, so temperature advection is occurring
    temp_global = xr.open_dataset(temp_path).sel(time=slice(str(y1), str(y2)))
    io_sst_weighted_smooth = (
        io_weighted_smooth  # use the smoothed version from SST processing
    )
    io_sst_smoothed_norm = io_sst_weighted_smooth.copy()
    io_sst_smoothed_norm.values = zscore(io_sst_weighted_smooth)
    io_sm_pos_phase = io_sst_smoothed_norm.where(io_sst_smoothed_norm > 0)
    temp_io_phase = temp_global.copy()
    temp_io_phase["time"] = io_sm_pos_phase.time
    temp_IO_pos = temp_global.where(io_sm_pos_phase).t2m.mean("time")
    temp_IO_pos_zonal_anom = temp_IO_pos - temp_IO_pos.mean("longitude")

    return (
        temp_reg_tropatl_masked,
        temp_reg_io_masked,
        malavg_temp,
        temp_IO_pos_zonal_anom,
    )


def process_soil_texture_data(soil_texture_path):
    """process GLDAS soil texture data for supplementary analysis"""
    print("Processing soil texture data...")

    # GLDAS
    gldas = xr.open_dataset(soil_texture_path).drop("time_bnds").mean("time")

    # for some reason the values are floats (~1e-5 away from int). Round to the int values which have physical meaning mapping 2 the categories
    gldas["GLDAS_soiltex"].values = np.round(gldas["GLDAS_soiltex"].values).astype(int)
    gldas = gldas.where(gldas.GLDAS_soiltex != 0, drop=True)

    # the GLDAS data uses "pixel-centered" coordinates while the clipping
    # operation expects "grid-centered" coordinates, causing apparent rightward shift
    # convert from pixel-centered to grid-centered by shifting coordinates
    # by half a pixel (0.125 deg) to the left and down
    new_lons = np.arange(-180, 180, 0.25)
    new_lats = np.arange(-60, 90.25, 0.25)

    gldas_resampled = gldas.interp(lon=new_lons, lat=new_lats, method="nearest")

    # colors for map based. colors from crameri batlow https://doi.org/10.5281/zenodo.1243862
    colors_mapped_soiltex = {
        3: "#fdb9c2",
        6: "#f6a077",
        7: "#c49138",
        8: "#818232",
        9: "#48714f",
        10: "#1e5c62",
        12: "#0f3c5f",
    }

    # soil labels from GLDAS
    soil_labels = {
        1: "Sand",
        2: "Loamy Sand",
        3: "Sandy Loam",
        4: "Silt Loam",
        5: "Silt",
        6: "Loam",
        7: "Sandy Clay Loam",
        8: "Silty Clay Loam",
        9: "Clay Loam",
        10: "Sandy Clay",
        11: "Silty Clay",
        12: "Clay",
        13: "Organic Materials",
        14: "Water",
        15: "Bedrock",
        16: "Other",
    }

    values_soiltex = sorted(colors_mapped_soiltex.keys())
    colors_soiltex = [colors_mapped_soiltex[v] for v in values_soiltex]
    cmap_soiltex = ListedColormap(colors_soiltex)
    boundaries_soiltex = [v - 0.5 for v in values_soiltex] + [values_soiltex[-1] + 0.5]
    norm_soiltex = BoundaryNorm(boundaries_soiltex, len(colors_soiltex))

    return gldas_resampled, colors_mapped_soiltex, soil_labels


def save_preprocessed_data(output_dir, **data_dict):
    """save all preprocessed data to appropriate file formats"""
    print(f"Saving preprocessed data to {output_dir}")

    for name, data in data_dict.items():
        try:
            if isinstance(data, xr.DataArray) or isinstance(data, xr.Dataset):
                # save xarray objects as netcdf
                filepath = os.path.join(output_dir, f"{name}.nc")
                data.to_netcdf(filepath)
                print(f"saved {name}.nc")

            elif isinstance(data, (pd.DataFrame, gpd.GeoDataFrame)):
                # save dataframes as csv (or parquet for geodataframes to preserve geometry)
                if isinstance(data, gpd.GeoDataFrame):
                    filepath = os.path.join(output_dir, f"{name}.parquet")
                    data.to_parquet(filepath)
                    print(f"saved {name}.parquet")
                else:
                    filepath = os.path.join(output_dir, f"{name}.csv")
                    data.to_csv(filepath, index=False)
                    print(f"saved {name}.csv")

            elif isinstance(data, np.ndarray):
                # save numpy arrays as npy
                filepath = os.path.join(output_dir, f"{name}.npy")
                np.save(filepath, data)
                print(f"saved {name}.npy")

            elif isinstance(data, dict):
                # save dictionaries as json where possible
                import json

                filepath = os.path.join(output_dir, f"{name}.json")
                with open(filepath, "w") as f:
                    json.dump(data, f, indent=2, default=str)
                print(f"saved {name}.json")

            else:
                # fallback to pickle for other data types
                import pickle

                filepath = os.path.join(output_dir, f"{name}.pkl")
                with open(filepath, "wb") as f:
                    pickle.dump(data, f)
                print(f"saved {name}.pkl")

        except Exception as e:
            print(f"failed to save {name}: {str(e)}")

    print("all preprocessed data saved successfully!")


def main():
    """main preprocessing pipeline"""
    print("Starting malaria-climate data preprocessing pipeline...")

    # create output directory
    output_dir = create_output_directory(OUTPUT_DIR)

    # load malaria data
    mal, y1, y2, casevar = load_malaria_data(MALARIA_DATA_PATH)

    # load supporting data
    malawi_shape, dist_shape, lake_malawi, etopo01 = load_supporting_data(
        MALAWI_SHAPE_PATH, DISTRICT_SHAPE_PATH, LAKE_MALAWI_PATH, ETOPO_PATH
    )

    # compute malaria statistics
    (
        mal,
        malavg,
        malavgscale,
        nattot,
        cases_pivot,
        pc1,
        pc2,
        sm_pc1_cases_24,
        sm_pc2_cases_24,
        explained_variance,
    ) = compute_malaria_statistics(mal, dist_shape, casevar)

    # process SST data
    (
        tropatl_weighted_nosmooth,
        io_sst_weighted_nosmooth,
        tropatl_anom_nosmooth,
        io_anom_nosmooth,
        tropatl_weighted_smooth,
        io_weighted_smooth,
        correlation_map_pc1,
        correlation_map_pc2,
    ) = process_sst_data(SST_PATH, y1, y2, pc1, pc2)

    # process wind data
    (u, v, u_mean, v_mean, u_reg_tropatl, v_reg_tropatl, u_reg_io, v_reg_io) = (
        process_wind_data(WIND_PATH, y1, y2, tropatl_anom_nosmooth, io_anom_nosmooth)
    )

    # process geopotential data
    (
        geop_850_height_mean,
        reg_geop_300_anom_tropatl,
        reg_geop_300_anom_io,
        reg_geop_850_anom_tropatl,
        reg_geop_850_anom_io,
    ) = process_geopotential_data(
        GEOP_300_PATH,
        GEOP_850_PATH,
        y1,
        y2,
        tropatl_anom_nosmooth,
        io_anom_nosmooth,
        tropatl_weighted_smooth,
        tropatl_weighted_nosmooth,
    )

    # process precipitation data
    (precip_reg_masked_TA, precip_reg_masked_IO, malavg_precip) = (
        process_precipitation_data(
            PRECIP_PATH,
            y1,
            y2,
            malawi_shape,
            tropatl_weighted_smooth,
            io_weighted_smooth,
            dist_shape,
            mal,
        )
    )

    # process VIMFD data
    reg_vimdf_TA_smoothed, reg_vimdf_IO_smoothed = process_vimfd_data(
        VIMFD_PATH, y1, y2, tropatl_anom_nosmooth, io_anom_nosmooth
    )

    # process soil moisture data
    (
        soil_reg_tropatl_masked,
        soil_reg_io_masked,
        correlation_map_soilm_pc1,
        correlation_map_soilm_pc2,
        malavg_soilm,
        sm_soilm_pc1,
        sm_soilm_pc2,
    ) = process_soil_moisture_data(
        GRACE_PATH,
        GRACE_MALAWI_PATH,
        y1,
        y2,
        tropatl_anom_nosmooth,
        io_anom_nosmooth,
        sm_pc1_cases_24,
        sm_pc2_cases_24,
        casevar,
        dist_shape,
        mal,
    )

    # process temperature data
    (
        temp_reg_tropatl_masked,
        temp_reg_io_masked,
        malavg_temp,
        temp_IO_pos_zonal_anom,
    ) = process_temperature_data(
        TEMP_PATH,
        y1,
        y2,
        malawi_shape,
        tropatl_anom_nosmooth,
        io_anom_nosmooth,
        dist_shape,
        mal,
        io_weighted_smooth,
    )

    # process soil texture data
    gldas_resampled, colors_mapped_soiltex, soil_labels = process_soil_texture_data(
        SOIL_TEXTURE_PATH
    )

    # prepare all data for savingg
    data_to_save = {
        # basic parameters
        "y1": y1,
        "y2": y2,
        "casevar": casevar,
        # geographic
        "malawi_shape": malawi_shape,
        "dist_shape": dist_shape,
        "lake_malawi": lake_malawi,
        "etopo01": etopo01,
        # malaria
        "mal": mal,
        "malavg": malavg,
        "malavgscale": malavgscale,
        "nattot": nattot,
        "cases_pivot": cases_pivot,
        "pc1": pc1,
        "pc2": pc2,
        "sm_pc1_cases_24": sm_pc1_cases_24,
        "sm_pc2_cases_24": sm_pc2_cases_24,
        "explained_variance": explained_variance,
        # SST
        "tropatl_weighted_nosmooth": tropatl_weighted_nosmooth,
        "io_sst_weighted_nosmooth": io_sst_weighted_nosmooth,
        "tropatl_anom_nosmooth": tropatl_anom_nosmooth,
        "io_anom_nosmooth": io_anom_nosmooth,
        "tropatl_weighted_smooth": tropatl_weighted_smooth,
        "io_weighted_smooth": io_weighted_smooth,
        "correlation_map_pc1": correlation_map_pc1,
        "correlation_map_pc2": correlation_map_pc2,
        # wind
        "u": u,
        "v": v,
        "u_mean": u_mean,
        "v_mean": v_mean,
        "u_reg_tropatl": u_reg_tropatl,
        "v_reg_tropatl": v_reg_tropatl,
        "u_reg_io": u_reg_io,
        "v_reg_io": v_reg_io,
        # geopotential
        "geop_850_height_mean": geop_850_height_mean,
        "reg_geop_300_anom_tropatl": reg_geop_300_anom_tropatl,
        "reg_geop_300_anom_io": reg_geop_300_anom_io,
        "reg_geop_850_anom_tropatl": reg_geop_850_anom_tropatl,
        "reg_geop_850_anom_io": reg_geop_850_anom_io,
        # precipitation
        "precip_reg_masked_TA": precip_reg_masked_TA,
        "precip_reg_masked_IO": precip_reg_masked_IO,
        "malavg_precip": malavg_precip,
        # VIMFD
        "reg_vimdf_TA_smoothed": reg_vimdf_TA_smoothed,
        "reg_vimdf_IO_smoothed": reg_vimdf_IO_smoothed,
        # soil moisture
        "soil_reg_tropatl_masked": soil_reg_tropatl_masked,
        "soil_reg_io_masked": soil_reg_io_masked,
        "correlation_map_soilm_pc1": correlation_map_soilm_pc1,
        "correlation_map_soilm_pc2": correlation_map_soilm_pc2,
        "malavg_soilm": malavg_soilm,
        "sm_soilm_pc1": sm_soilm_pc1,
        "sm_soilm_pc2": sm_soilm_pc2,
        # temperature
        "temp_reg_tropatl_masked": temp_reg_tropatl_masked,
        "temp_reg_io_masked": temp_reg_io_masked,
        "malavg_temp": malavg_temp,
        "temp_IO_pos_zonal_anom": temp_IO_pos_zonal_anom,
        # soil texture
        "gldas_resampled": gldas_resampled,
        "colors_mapped_soiltex": colors_mapped_soiltex,
        "soil_labels": soil_labels,
    }

    # save all preprocessed data
    save_preprocessed_data(output_dir, **data_to_save)

    print("preprocessing pipeline completed successfully!")


if __name__ == "__main__":
    main()
