"""
Utilities for modeling and analyzing snow cover trends
Rainey Aberle
2023
"""

import math
import glob
import os
from tqdm.auto import tqdm
import numpy as np
import xarray as xr
import rioxarray as rxr
import rasterio


def convert_wgs_to_utm(lon: float, lat: float):
    """
    Return best UTM epsg-code based on WGS84 lat and lon coordinate pair

    Parameters
    ----------
    lon: float
        longitude coordinate
    lat: float
        latitude coordinate

    Returns
    ----------
    epsg_code: str
        optimal UTM zone, e.g. "32606"
    """
    utm_band = str(int((np.floor((lon + 180) / 6) % 60) + 1))
    if len(utm_band) == 1:
        utm_band = '0' + utm_band
    if lat >= 0:
        epsg_code = 'EPSG:326' + utm_band
        return epsg_code
    epsg_code = 'EPSG:327' + utm_band
    return epsg_code


def determine_subregion_name_color(o1, o2):
    if (o1 == 1.0) and (o2 == 1.0):
        subregion_name, color = 'Brooks Range', 'c'
    elif (o1 == 1.0) and (o2 == 2.0):
        subregion_name, color = 'Alaska Range', '#1f78b4'
    elif (o1 == 1.0) and (o2 == 3.0):
        subregion_name, color = 'Aleutians', '#6d9c43'
    elif (o1 == 1.0) and (o2 == 4.0):
        subregion_name, color = 'W. Chugach Mtns.', '#264708'
    elif (o1 == 1.0) and (o2 == 5.0):
        subregion_name, color = 'St. Elias Mtns.', '#fb9a99'
    elif (o1 == 1.0) and (o2 == 6.0):
        subregion_name, color = 'N. Coast Ranges', '#e31a1c'
    elif (o1 == 2.0) and (o2 == 1.0):
        subregion_name, color = 'N. Rockies', '#cab2d6'
    elif (o1 == 2.0) and (o2 == 2.0):
        subregion_name, color = 'N. Cascades', '#fdbf6f'
    elif (o1 == 2.0) and (o2 == 3.0):
        subregion_name, color = 'C. Rockies', '#9657d9'
    elif (o1 == 2.0) and (o2 == 4.0):
        subregion_name, color = 'S. Cascades', '#ff7f00'
    elif (o1 == 2.0) and (o2 == 5.0):
        subregion_name, color = 'S. Rockies', '#6a3d9a'
    else:
        subregion_name = 'O1:' + o1 + ' O2:' + o2
        color = 'k'

    return subregion_name, color


def load_snow_cover_stats(fn):
    ds = xr.open_zarr(fn)
    cols = ['AAR', 'glacier_area', 'ice_area', 'SLA', 'SLA_lower_bound', 'SLA_upper_bound', 
            'snow_area', 'water_area']
    ds = ds[cols] # don't need classification band (much larger)
    # convert dask arrays to numpy arrays for easier computations
    for col in cols:
        ds[col] = ds[col].compute()
    ds['SLA'] = xr.where(ds['SLA']==0, np.nan, ds['SLA'])
    return ds


def calculate_image_snow_cover_stats(fn, epsg_utm, reference_raster, dem, dem_min, dem_max, aoi_area):
    # Load classified image
    ds = xr.open_dataset(fn)
    ds = ds.rio.write_crs("EPSG:4326")
    ds = ds.rio.reproject_match(reference_raster)
    dx = np.nanmean(ds.x.values[1:] - ds.x.values[0:-1])

    # Calculate land cover areas
    snow_area = len(np.ravel(ds.classified.data[ds.classified.data <= 2])) * (dx ** 2)
    ice_area = len(np.ravel(ds.classified.data[ds.classified.data == 3])) * (dx ** 2)
    rock_area = len(np.ravel(ds.classified.data[ds.classified.data == 4])) * (dx ** 2)
    water_area = len(np.ravel(ds.classified.data[ds.classified.data == 5])) * (dx ** 2)
    glacier_area = ice_area + snow_area

    # Calculate AAR
    if glacier_area==0:
        aar, sla, sla_lower_bound, sla_upper_bound = np.nan, np.nan, np.nan, np.nan
    else:
        aar = snow_area / glacier_area
        if aar==1:
            sla, sla_lower_bound, sla_upper_bound = dem_min, dem_min, dem_min
        elif aar == 0:
            sla, sla_lower_bound, sla_upper_bound = dem_max, dem_max, dem_max
        else:

            # Calculate snowline altitude
            dem_adj = dem.rio.reproject_match(ds)
            elevations = dem_adj.data.ravel()
            sla = np.nanquantile(elevations, 1 - aar)

            # Calculate snowline altitude upper and lower bounds
            # lower bound
            snow_below_sla = xr.where((dem_adj < sla) & (ds.classified <= 2), 1, 0).data.ravel().sum() * (dx**2)
            sla_lower_percentile = aar - (snow_below_sla / aoi_area)
            sla_lower_percentile = np.clip(sla_lower_percentile, 0, 1)
            sla_lower_bound = np.nanquantile(elevations, sla_lower_percentile)
            # upper bound
            snow_free_above_sla = xr.where((dem_adj > sla) & (ds.classified > 2), 1, 0).data.ravel().sum() * (dx**2)
            sla_upper_percentile = aar + (snow_free_above_sla / aoi_area)
            sla_upper_percentile = np.clip(sla_upper_percentile, 0, 1)
            sla_upper_bound = np.nanquantile(elevations, sla_upper_percentile)
    
    # Reproject back to lat/lon
    ds = ds.rio.write_crs(epsg_utm)
    ds = ds.rio.reproject("EPSG:4326")

    # Add stats as data variables
    for name, area in [['snow_area', snow_area],
                    ['ice_area', ice_area],
                    ['rock_area', rock_area],
                    ['water_area', water_area],
                    ['glacier_area', glacier_area]]:
        ds[name] = (('time'), [area])
        ds[name].attrs['units'] = 'meters squared'
    ds['AAR'] = (('time'), [aar])
    ds['SLA'] = (('time'), [sla])
    ds['SLA_lower_bound'] = (('time'), [sla_lower_bound])
    ds['SLA_upper_bound'] = (('time'), [sla_upper_bound])

    return ds


def save_classified_image(ds, rgi_id, epsg_utm, out_fn):
    # Convert classification band to int data format
    if 'classified' in ds.data_vars:
        ds = ds.rename_vars({'classified': 'classification'})
    ds['classification'] = ds['classification'].fillna(0).astype(np.uint8)

    # Assign global attributes
    ds.attrs['title'] = 'Land cover classifications and snow cover stastics'
    ds.attrs['description'] = 'Land cover classifications and snow cover statistics constructed using the glacier-snow-cover-mapping package (https://github.com/RaineyAbe/glacier-snow-cover-mapping).'
    ds.attrs['references'] = 'doi:10.1029/2025GL115523'
    ds.attrs['horizontal_CRS'] = 'WGS84 (EPSG:4326)'
    ds.attrs['vertical_CRS'] = 'EGM96 geoid (EPSG:5773)'
    ds.attrs['date_modified'] = '2025-06-07'
    ds.attrs['time_coverage_start'] = '2013-05-01'
    ds.attrs['time_coverage_end'] = '2023-10-31'
    ds.attrs['spatial_resolution'] = '10-30m, depending on the image source used to construct the classified image.'
    ds.attrs['RGIId'] = rgi_id

    # Assign variable-specific attributes
    ds['classification'].attrs['long_name'] = 'land cover classification'
    ds['classification'].attrs['classes'] = '1 = snow, 2 = shadowed snow, 3 = ice, 4 = rock/debris, 5 = water'
    ds['classification'].attrs['NoData'] = 0
    ds['AAR'].attrs['long_name'] = 'accumulation area ratio'
    ds['AAR'].attrs['units'] = 'unitless'
    ds['SLA'].attrs['long_name'] = 'snowline altitude'
    ds['SLA'].attrs['units'] = 'meters above sea level'
    ds['SLA_lower_bound'].attrs['long_name'] = 'snowline altitude lower bound'
    ds['SLA_lower_bound'].attrs['units'] = 'meters above sea level'
    ds['SLA_upper_bound'].attrs['long_name'] = 'snowline altitude upper bound'
    ds['SLA_upper_bound'].attrs['units'] = 'meters above sea level'
    ds['glacier_area'].attrs['units'] = 'meters above sea level'
    ds['glacier_area'].attrs['long_name'] = 'glacier area'
    ds['glacier_area'].attrs['description'] = 'snow area + ice area'
    for dv in ['snow_area', 'ice_area', 'rock_area', 'water_area']:
        ds[dv].attrs['units'] = 'meters squared'
        ds[dv].attrs['long_name'] = (dv).split('_')[0] + ' cover area'

    # Sort data variables alphabetically
    ds = ds[sorted(ds.data_vars.keys())]

    # Make sure no empty "spatial_ref" coordinate
    if 'spatial_ref' in ds.coords:
        ds = ds.drop_vars('spatial_ref')

    # Save to file
    if '.nc' in out_fn:
        ds.to_netcdf(out_fn)
        print('Compiled dataset saved to file:', out_fn)
    elif '.zarr' in out_fn:
        ds.to_zarr(out_fn)
        print('Compiled dataset saved to file:', out_fn)


def compile_classified_image_files(scm_path, rgi_id, aoi):
    # Define output file name
    out_fn = os.path.join(scm_path, 'study-sites', rgi_id, f"{rgi_id}_classifications.zarr")
    if os.path.exists(out_fn):
        return
    
    # Get classified image file names
    classified_fns = sorted(glob.glob(os.path.join(scm_path, 'study-sites', rgi_id, 'classified', f"*classified.nc")))

    # get optimal UTM zone
    cen_lon, cen_lat = aoi.geometry[0].centroid.coords.xy
    cen_lon, cen_lat = cen_lon[0], cen_lat[0]
    epsg_utm = convert_wgs_to_utm(cen_lon, cen_lat)
    aoi_utm = aoi.to_crs(epsg_utm)
    aoi_area = aoi_utm.area.values[0]

    # Create a 10 m grid for all image files
    xmin, ymin, xmax, ymax = aoi_utm.total_bounds
    resolution = 10
    cols = int((xmax - xmin) / resolution)
    rows = int((ymax - ymin) / resolution)
    transform = rasterio.transform.from_bounds(xmin, ymin, xmax, ymax, cols, rows)
    x = np.arange(xmin, xmax, resolution)
    y = np.arange(ymax, ymin, -resolution)
    reference_raster = xr.DataArray(
        np.nan * np.zeros((len(y), len(x))),
        coords={"y": y, "x": x},
        dims=("y", "x"),
        attrs={"transform": transform, "crs": epsg_utm}
    )

    # Load DEM
    dem_fn = sorted(glob.glob(os.path.join(scm_path, 'study-sites', rgi_id, 'DEMs', '*.tif')))[0]
    dem = rxr.open_rasterio(dem_fn)
    dem = dem.isel(band=0)
    dem = dem.rio.reproject_match(reference_raster)
    dem = dem.rio.clip(aoi_utm.geometry)
    dem = xr.where((dem > 1e4) | (dem < -1e3), np.nan, dem)
    dem = dem.rio.write_crs(epsg_utm)
    dem_min = float(dem.min().data)
    dem_max = float(dem.max().data)

    # If glacier area > 1000 km2, must be saved in chunks then compiled
    if aoi['Area'].values[0] < 1e3:

        # Initialize list of datasets
        ds_list = []

        # Iterate over files
        for j, fn in enumerate(tqdm(classified_fns)):
            # process file
            ds = calculate_image_snow_cover_stats(fn, epsg_utm, reference_raster, dem, dem_min, dem_max, aoi_area)
            ds_list.append(ds)

        # Concatenate into single dataset
        ds_full = xr.concat(ds_list, dim='time', coords='different')

        # Save to file
        save_classified_image(ds_full, rgi_id, epsg_utm, out_fn)
    
    else:
        # split classified images into 10 sets
        idx = np.linspace(0, len(classified_fns), num=10).astype(int)
        for i in np.arange(len(idx)-1):
            # Define output file name
            out_fn = os.path.join(scm_path, 'study-sites', rgi_id, f"{rgi_id}_classifications_{i}.zarr")
            if os.path.exists(out_fn):
                continue

            # Initialize list of datasets
            ds_list = []
            classified_fns_idx = classified_fns[idx[i]:idx[i+1]]

            # Iterate over files
            for j, fn in enumerate(tqdm(classified_fns_idx)):
                # process file
                ds = calculate_image_snow_cover_stats(fn, epsg_utm, reference_raster, dem, dem_min, dem_max, aoi_area)
                ds_list.append(ds)

            # Concatenate into single dataset
            ds_full = xr.concat(ds_list, dim='time', coords='different')

            # Save to file
            save_classified_image(ds_full, rgi_id, epsg_utm, out_fn)

        # combine all!
        out_fns = glob.glob(os.path.join(scm_path, 'study-sites', rgi_id, f"{rgi_id}_classifications_*.zarr"))
        ds_list = []
        for out_fn in out_fns:
            ds = xr.open_zarr(out_fn)
            ds_list.append(ds)
        ds_full = xr.concat(ds_list, dim='time', coords='different')
        ds_full = ds_full.sortby('time')
        ds_full = ds_full.chunk({'time': len(ds_full.time.data)})
        final_out_fn = os.path.join(scm_path, 'study-sites', rgi_id, f"{rgi_id}_classifications.zarr")
        save_classified_image(ds_full, rgi_id, epsg_utm, final_out_fn)
    

def calculate_hypsometric_index(dem_fn, aoi):
    """
    Calculate Hypsometric Index using an input DEM file name and area of interest shapefile. 
    Based on Jiskoot et al. (2009): https://doi.org/10.3189/172756410790595796

    Parameters
    ----------
    dem_fn: str | Path
        file name of the DEM
    aoi: geopandas.GeoDataFrame
        area of interest (AOI) with "geometry" and "CRS" set

    Returns
    ----------
    hi: float
        hypsometric index
    hi_category: str
        hypsometric index category ("Very bottom heavy", "Bottom heavy", "Equidimensional", 
        "Top heavy", or "Very top heavy)
    """
    
    # load DEM as DataArray
    dem = rxr.open_rasterio(dem_fn)
    # reproject DEM to AOI CRS
    dem = dem.rio.reproject('EPSG:'+str(aoi.crs.to_epsg()))
    # clip DEM to AOI
    try:
        dem_aoi = dem.rio.clip(aoi.geometry, aoi.crs)
    except:
        return 'N/A', 'N/A'
    # convert to dataset
    dem_aoi_ds = dem_aoi.to_dataset(name='elevation')
    # set no data values to NaN
    dem_aoi_ds = xr.where((dem_aoi_ds > 1e38) | (dem_aoi_ds <= -9999), np.nan, dem_aoi_ds)
    # check that there is data after removing no data values
    if np.isnan(dem_aoi_ds.elevation.data).all():
        return 'N/A', 'N/A'
    # calculate elevation statistics
    h_max = np.nanmax(np.ravel(dem_aoi_ds.elevation.data))
    h_min = np.nanmin(np.ravel(dem_aoi_ds.elevation.data))
    h_med = np.nanmedian(np.ravel(dem_aoi_ds.elevation.data))
    # calculate HI, where HI = (H_max - H_med) / (H_med - H_min). If 0 < HI < 1, HI = -1/HI.
    hi = (h_max - h_med) / (h_med - h_min)
    if (0 < hi) and (hi < 1):
        hi = -1 / hi
    # determine HI category
    if hi <= -1.5:
        hi_category = 'Very top heavy'
    elif (hi > -1.5) and (hi <= -1.2):
        hi_category = 'Top heavy'
    elif (hi > -1.2) and (hi <= 1.2):
        hi_category = 'Equidimensional'
    elif (hi > 1.2) and (hi <= 1.5):
        hi_category = 'Bottom heavy'
    elif hi > 1.5:
        hi_category = 'Very bottom heavy'

    return hi, hi_category


def reduce_df_memory_usage(df, verbose=True):
    """
    Reduce memory usage in pandas.DataFrame
    From Bex T (2021): https://towardsdatascience.com/6-pandas-mistakes-that-silently-tell-you-are-a-rookie-b566a252e60d

    Parameters
    ----------
    df: pandas.DataFrame
        input dataframe
    verbose: bool
        whether to output verbage (default=True)

    Returns
    ----------
    df: pandas.DataFrame
        output dataframe with reduced memory usage
    """
    numerics = ["int8", "int16", "int32", "int64", "float16", "float32", "float64"]
    start_mem = df.memory_usage().sum() / 1024 ** 2
    for col in df.columns:
        col_type = df[col].dtypes
        if col_type in numerics:
            c_min = df[col].min()
            c_max = df[col].max()
            if str(col_type)[:3] == "int":
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)
            else:
                #                if (
                #                    c_min > np.finfo(np.float16).min
                #                    and c_max < np.finfo(np.float16).max
                #                ):
                #                    df[col] = df[col].astype(np.float16) # float16 not compatible with linalg
                if (  # elif (
                        c_min > np.finfo(np.float32).min
                        and c_max < np.finfo(np.float32).max
                ):
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
    end_mem = df.memory_usage().sum() / 1024 ** 2
    if verbose:
        print(
            "pandas.DataFrame memory usage decreased to {:.2f} Mb ({:.1f}% reduction)".format(
                end_mem, 100 * (start_mem - end_mem) / start_mem
            )
        )
    return df
