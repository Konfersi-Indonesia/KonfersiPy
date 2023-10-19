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
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
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

def reproject_shp_gdal(infile, outfile, targetprj):
    driver = ogr.GetDriverByName("ESRI Shapefile") 
    dataSource = driver.Open(infile, 1)
    layer = dataSource.GetLayer()
    sourceprj = layer.GetSpatialRef()
    transform = osr.CoordinateTransformation(sourceprj, targetprj)
    outDriver = ogr.GetDriverByName("Esri Shapefile")
    outDataSource = outDriver.CreateDataSource(outfile)
    outlayer = outDataSource.CreateLayer('', targetprj, ogr.wkbPolygon)
    outlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    i = 0
    for feature in layer:
        transformed = feature.GetGeometryRef()
        transformed.Transform(transform)
        geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
        defn = outlayer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', i)
        feat.SetGeometry(geom)
        outlayer.CreateFeature(feat) 
        i += 1
        feat = None
    
    print("Driver Reproyeksi: ", driver)
    print("Data Source: ", dataSource)
    print("Layer: ", layer)
    print("Source Projection: ", sourceprj)
    print("Transform: ", transform)
    print("Driver Hasil: ", outDriver)
    print("Hasil Data: ", outDataSource)
    print("Hasil Layer: ", outlayer)
    print("Feature: ", feature)

def array2raster(array, geoTransform, projection, filename):
    pixels_x = array.shape[1]
    pixels_y = array.shape[0]
    
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(
        filename,
        pixels_x,
        pixels_y,
        1,
        gdal.GDT_Float64, )
    dataset.SetGeoTransform(geoTransform)
    dataset.SetProjection(projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()
    
    print("Driver: ", driver)
    print("Dataset: ", dataset)
    return dataset, dataset.GetRasterBand(1)

def clip_raster(filename, shp):
    inraster = rasterio.open(filename)
    extent_geojson = mapping(shp['geometry'][0])
    clipped, crop_affine = mask(inraster, 
                                shapes=[extent_geojson], 
                                nodata = np.nan,
                                crop=True)
    clipped_meta = inraster.meta.copy()
    clipped_meta.update({"driver": "GTiff",
                 "height": clipped.shape[0],
                 "width": clipped.shape[1],
                 "transform": crop_affine})
    cr_ext = rasterio.transform.array_bounds(clipped_meta['height'], 
                                            clipped_meta['width'], 
                                            clipped_meta['transform'])
    gt = crop_affine.to_gdal()
    
    print("Extent Geojson: ", extent_geojson)
    print("Clip: ", clipped)
    print("Crop Affine: ", crop_affine)
    print("Clipped Meta: ", clipped_meta)
    print("Crop Extent: ", cr_ext)
    print("Geotransform: ", gt)
    
    return clipped, clipped_meta, cr_ext, gt

def read_band_image(band, path):
    from osgeo import gdal
    a = path+band
    img = gdal.Open(glob.glob(a)[0])
    data = np.array(img.GetRasterBand(1).ReadAsArray())
    spatialRef = img.GetProjection()
    geoTransform = img.GetGeoTransform()
    targetprj = osr.SpatialReference(wkt = img.GetProjection())
    
    print("Data Band: ", data)
    print("Spasial Referensi: ", spatialRef)
    print("Geotransform: ", geoTransform)
    print("Target Proyeksi: ", targetprj)
    
    return data, spatialRef, geoTransform, targetprj

def reproject_shp_gdal(infile, outfile, targetprj):
    driver = ogr.GetDriverByName("ESRI Shapefile") 
    dataSource = driver.Open(infile, 1)
    layer = dataSource.GetLayer()
    sourceprj = layer.GetSpatialRef()
    transform = osr.CoordinateTransformation(sourceprj, targetprj)
    outDriver = ogr.GetDriverByName("Esri Shapefile")
    outDataSource = outDriver.CreateDataSource(outfile)
    outlayer = outDataSource.CreateLayer('', targetprj, ogr.wkbPolygon)
    outlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    i = 0
    for feature in layer:
        transformed = feature.GetGeometryRef()
        transformed.Transform(transform)
        geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
        defn = outlayer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', i)
        feat.SetGeometry(geom)
        outlayer.CreateFeature(feat) 
        i += 1
        feat = None
    
    print("Driver Reproyeksi: ", driver)
    print("Data Source: ", dataSource)
    print("Layer: ", layer)
    print("Source Projection: ", sourceprj)
    print("Transform: ", transform)
    print("Driver Hasil: ", outDriver)
    print("Hasil Data: ", outDataSource)
    print("Hasil Layer: ", outlayer)
    print("Feature: ", feature)

def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)

def timeIndexToDatetime(baseTime,times):
    newTimes=[]
    for ts in times:
        newTimes.append(baseTime+datetime.timedelta(seconds=ts))

    return newTimes
    
def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and 
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray
