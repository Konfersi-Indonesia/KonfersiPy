from netCDF4 import Dataset
import netCDF4 as netcdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import xarray as xr
import pandas as pd
import glob
import os
import sys
import rasterio
from rasterio import plot
from rasterio.enums import Resampling
from rasterio.plot import show
from rasterio.mask import mask
from osgeo import ogr, gdal, osr
import fiona
import shapely
from shapely.geometry import mapping
import geopandas as gpd
import math
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from mpl_toolkits import basemap
import matplotlib as mpl

def wind_stress(u, v, rho_air=1.22, cd=None):
    """Convert wind speed (u,v) to wind stress (Tx,Ty).

    It uses either wind-dependent or constant drag.

    Args:
        u, v: Wind vector components (m/s), 2d or 3d (for time series).
        rho_air: Density of air (1.22 kg/m^3).
        cd: Non-dimensional drag (wind-speed dependent).
            For constant drag use cd=1.5e-3.
    Notes:
        Function to compute wind stress from wind field data is based on Gill,
        (1982)[1]. Formula and a non-linear drag coefficient (cd) based on
        Large and Pond (1981)[2], modified for low wind speeds (Trenberth et
        al., 1990)[3]

        [1] A.E. Gill, 1982, Atmosphere-Ocean Dynamics, Academy Press.
        [2] W.G. Large & S. Pond., 1981,Open Ocean Measurements in Moderate
        to Strong Winds, J. Physical Oceanography, v11, p324-336.
        [3] K.E. Trenberth, W.G. Large & J.G. Olson, 1990, The Mean Annual
        Cycle in Global Ocean Wind Stress, J. Physical Oceanography, v20,
        p1742-1760.
    """
    w = np.sqrt(u**2 + v**2) # wind speed (m/s) 

    if not cd:
        # wind-dependent drag
        cond1 = (w<=1)
        cond2 = (w>1) & (w<=3)
        cond3 = (w>3) & (w<10)
        cond4 = (w>=10)
        cd = np.zeros_like(w)
        cd[cond1] = 2.18e-3 
        cd[cond2] = (0.62 + 1.56/w[cond2]) * 1e-3
        cd[cond3] = 1.14e-3
        cd[cond4] = (0.49 + 0.065*w[cond4]) * 1e-3

    Tx = rho_air * cd * w * u # zonal wind stress (N/m^2)
    Ty = rho_air * cd * w * v # meridional wind stress (N/m^2)
    return [Tx, Ty]


def wind_curl(u, v, x, y, ydim=0, xdim=1, tdim=2):
    """Calculate the curl of the wind vector. 
    
    Args:
        u, v: Wind vector components (m/s), 2d or 3d (for time series).
        x, y: Coordinates in lon/lat (degrees), 2d.

    Notes:
        Curl(u,v) = dv/dx - du/dy
        Units of frequency (1/s).
        The different constants come from oblateness of the ellipsoid.
    """
    dy = np.abs(y[1,0] - y[0,0]) # scalar in deg
    dx = np.abs(x[0,1] - x[0,0]) 
    dy *= 110575. # scalar in m
    dx *= 111303. * np.cos(y * np.pi/180) # array in m (varies w/lat)  # FIXME?
    # extend dimension for broadcasting (2d -> 3d)
    dvdx = []
    if u.ndim == 3:
        dx = np.expand_dims(dx, tdim)
        
    # grad[f(y,x), delta] = diff[f(y)]/delta, diff[f(x)]/delta 
    dudy = np.gradient(u, dy)[ydim] # (1/s)
    for i in range(0,len(v[:,0]),1):
        dvdx1 = np.gradient(v[i,:], dx[i,0])
        dvdx = np.concatenate((dvdx,dvdx1),axis=0)
    dvdx = np.reshape(dvdx,(len(v[:,0]),len(v[0,:])))
    curl = dvdx - dudy # (1/s)
    return curl


def wind_stress_curl(Tx, Ty, x, y, ydim=0, xdim=1, tdim=2):
    """Calculate the curl of wind stress (Tx, Ty).

    Args:
        Tx, Ty: Wind stress components (N/m^2), 2d or 3d (for time series)
        x, y: Coordinates in lon/lat (degrees), 2d.

    Notes:
        Curl(Tx,Ty) = dTy/dx - dTx/dy
        The different constants come from oblateness of the ellipsoid.
    """
    dy = np.abs(y[1,0] - y[0,0]) # scalar in deg
    dx = np.abs(x[0,1] - x[0,0]) 
    dy *= 110575. # scalar in m
    dx *= 111303. * np.cos(y * np.pi/180) # array in m (varies w/lat)
    # extend dimension for broadcasting (2d -> 3d)
    if Tx.ndim == 3:
        dx = np.expand_dims(dx, tdim)
    
    dTydx = []
    # grad[f(y,x), delta] = diff[f(y)]/delta, diff[f(x)]/delta 
    dTxdy = np.gradient(Tx, dy)[ydim] # (N/m^3)
    #dTydx = np.gradient(Ty, dx)[xdim] 
    
    for i in range(0,len(Ty[:,0]),1):
        dTydx1 = np.gradient(Ty[i,:], dx[i,0])
        dTydx = np.concatenate((dTydx,dTydx1),axis=0)
    dTydx = np.reshape(dTydx,(len(Ty[:,0]),len(Ty[0,:])))
    
    curl_tau = dTydx - dTxdy # (N/m^3)
    return curl_tau

def ekman_transport(tau_x, tau_y, y, rho_water=1028., tdim=2):
    """Calculate Ekman Mass Transport from Wind Curl.
    
    Args:
        tau_x = Wind Curl X, 2d or 3d (for time series)
        tau y = Wind Curl Y, 2d or 3d (for time series)
        y = Latitude grid, 2d
    """
    
    omega = 7.292115e-5 # rotation rate of the Earth (rad/s)
    f = np.abs(2 * omega * np.sin(y * np.pi/180)) # (rad/s)

    if f.shape != tau_y.shape:
        f = np.expand_dims(f, tdim)

    emtx = (tau_y)/(rho_water*f)
    emty = (tau_x)/(rho_water*f)

    return emtx,emty
def ekman_pumping(curl_tau, y, rho_water=1028., tdim=2):
    """Calculate Ekman pumping from wind-stress curl.

    Args:
        curl_tau: Wind stress curl (N/m^3), 2d or 3d (for time series).
        y: Latitude grid (degrees), 2d.

    Notes:
        We = Curl(tau)/rho*f (vertical velocity in m/s).
        f = Coriolis frequency (rad/s), latitude dependent.
        rho = Ocean water density (1028 kg/m^3).
    """
    # Coriolis frequency
    omega = 7.292115e-5 # rotation rate of the Earth (rad/s)
    f = np.abs(2 * omega * np.sin(y * np.pi/180)) # (rad/s)

    # Expand dimension for broadcasting (2d -> 3d)
    if f.shape != curl_tau.shape:
        f = np.expand_dims(f, tdim)

    # Ekman pumping
    We = curl_tau / (rho_water * f) # vertical velocity (m/s)
    return We
