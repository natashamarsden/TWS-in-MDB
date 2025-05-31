from shapely.geometry import mapping
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray
from scipy.stats import sem, t
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates


def mask_data(data, boundary):

    """Returns a clipped dataset of the original data.
    
    Args:
        data (xarray.DataArray): an array of values with corresponding spatial dimensions labelled "longitude" and "latitude".
        boundary (geopandas.geodataframe): a polygon object contained in a geodataframe.

    Returns:
        clipped_data (xarray.DataArray): a masked version of data which only contains entries with spatial coordinates inside boundary."""

    # set spatial dimensions to long and lat components
    data.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)

    # set coordinate system
    data.rio.write_crs("epsg:4326", inplace=True)

    # perform masking
    clipped_data = data.rio.clip(boundary.geometry.apply(mapping), boundary.crs, drop=False)

    # returns masked dataset
    return clipped_data

def calc_area_weights(latitudes, longitudes, radius=6371000):
    """
    Calculate area (in m²) of each lat/lon grid cell using the colatitude formula (as provided by Jae-Seung).
    """

    delta_lat = np.abs(np.diff(latitudes).mean())
    delta_lon = np.abs(np.diff(longitudes).mean())
        
    # Convert to radians
    dtheta = np.deg2rad(delta_lat)
    dphi = np.deg2rad(delta_lon)
    colatitudes = 90 - latitudes
    colat_rad = np.deg2rad(colatitudes)
    

    # compute area weights for each latitude
    area_weights_1d = radius ** 2 * np.sin(colat_rad) * dtheta * dphi # shape: (n_lat,)
    area_weights_2d = np.tile(area_weights_1d[:, np.newaxis], (1, len(longitudes)))  # shape: (n_lat, n_lon)

    # Convert to xarray DataArray for compatibility
    area_weights = xarray.DataArray(
        area_weights_2d,
        coords={"latitude": latitudes, "longitude": longitudes},
        dims=["latitude", "longitude"]
    )

    return area_weights

# Area-weighted average
def weighted_mean(data, weights):
    weighted_sum = (data * weights).sum(dim=["latitude", "longitude"], skipna=True)
    sum_of_weights = weights.where(~data.isnull()).sum(dim=["latitude", "longitude"])
    return weighted_sum / sum_of_weights

def process_data(model_data, boundary_shape):

    # data must be ET or precip data over period from 15th of prior month to 15th of this month, in m

    # mask data
    clipped_data = mask_data(data = model_data, boundary = boundary_shape)
    
    # Get lat/lon
    latitudes = clipped_data.latitude.values
    longitudes = clipped_data.longitude.values
    
    # Calculate area weights
    area_weights = calc_area_weights(latitudes, longitudes, radius=6371000)

    # get mean data over spatial boundary
    mean_data = weighted_mean(clipped_data, area_weights)
    
    # convert to numpy array
    mean_data = mean_data.values
    
    return mean_data

def process_GRACE(data_GRACE, boundary_shape):
    
    # mask data
    masked_data = mask_data(data=data_GRACE, boundary = boundary_shape)
    
    # Get lat/lon
    latitudes = masked_data.latitude.values
    longitudes = masked_data.longitude.values

    # Calculate area weights
    area_weights = calc_area_weights(latitudes, longitudes, radius=6371000)
    
    # Get mean over spatial boundary
    TWS = weighted_mean(masked_data, area_weights)/1000

    return TWS

def mask_data_modis(data, boundary):

    """Returns a clipped dataset of the original data.
    
    Args:
        data (xarray.DataArray): an array of values with corresponding spatial dimensions labelled "longitude" and "latitude".
        boundary (geopandas.geodataframe): a polygon object contained in a geodataframe.

    Returns:
        clipped_data (xarray.DataArray): a masked version of data which only contains entries with spatial coordinates inside boundary."""

    # set spatial dimensions to long and lat components
    data.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)

    # set coordinate system
    data.rio.write_crs("epsg:4326", inplace=True)

    # perform masking
    clipped_data = data.rio.clip(boundary.geometry.apply(mapping), boundary.crs, drop=True)

    # returns masked dataset
    return clipped_data

def calc_spatial_avg(model_data, boundary_shape):

    # mask model data to region
    clipped_data = mask_data(data = model_data, boundary = boundary_shape)

    # Get lat/lon
    latitudes = clipped_data.latitude.values
    longitudes = clipped_data.longitude.values

    # Calculate area weights
    area_weights = calc_area_weights(latitudes, longitudes, radius=6371000)

    # get mean data over spatial boundary
    mean_data = weighted_mean(clipped_data, area_weights)

    dates = clipped_data.time

    return mean_data, dates

def ensemble_graph(series_list, GRACE_data, dates, title, precip=None, xlim1=None, xlim2=None):
    
    min_len = min(len(ts) for ts in series_list)
    dates = dates[:min_len]
    GRACE_data = GRACE_data[:min_len]
    
    if precip is not None:
        precip = precip[:min_len]
    
    stacked = np.array([ts[:min_len] for ts in series_list])
    data = stacked.T

    mean_series = np.nanmean(data, axis=1)
    n = np.sum(~np.isnan(data), axis=1)
    stderr = sem(data, axis=1, nan_policy='omit')
    ci = t.ppf(0.975, df=n-1) * stderr  # 95% CI

    # Set up figure with two vertically stacked subplots
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)

    if precip is not None:

        # Convert dates to matplotlib internal float format
        dates_num = mdates.date2num(dates)
    
        # Calculate bar width (roughly 80% of average difference between dates)
        if len(dates_num) > 1:
            width = 0.8 * (dates_num[1] - dates_num[0])
        else:
            width = 1  # fallback if only one date

        ax0 = fig.add_subplot(gs[0])
        ax0.bar(dates, precip, width=width, color='blue', alpha=0.6)
        ax0.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax0.set_ylabel('AWRA P (m)')
        if xlim1 and xlim2:
            ax0.set_xlim(xlim1, xlim2)
        ax0.grid(True)
    else:
        ax0 = None

    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    ax1.plot(dates, mean_series, label='Mean of LSM Combinations', color='blue', linewidth=0.75)
    ax1.plot(dates, GRACE_data, label='GRACE Delta TWS', color='red', linewidth=0.75)
    ax1.fill_between(dates, mean_series - ci, mean_series + ci, color='blue', alpha=0.2, label='95% CI')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Delta TWS (m)')
    if xlim1 and xlim2:
        ax1.set_xlim(xlim1, xlim2)
    ax1.legend()
    ax1.grid(True)

    fig.suptitle(title, fontsize=16)
    return fig

    #plt.tight_layout()
    #plt.show()