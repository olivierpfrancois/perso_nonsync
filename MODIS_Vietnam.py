
""" """

import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

import os, re
from datetime import datetime, timedelta
from osgeo import gdal
import MODIS_gedata_toolbox as md
reload(md)
import multiprocessing as mp
import contextlib as context
# from xml.dom import NoDataAllowedErr
# from contextlib import contextmanager
import numpy as np

sys.path.append('/home/olivierp/OlivierGithub/QGIS-scripts')
import image_proc_functions as img_proc

###########################################################################
###########################################################################
# #PARAMETERS

# #DIRECTORIES parameters
# Working directory
dst = '/home/olivierp/jde_coffee/MODIS/collection6/Vietnam'
# Folder to use for temporary files (should be empty)
tempDir = '/home/olivierp/jde_coffee/Temp'
# Destination folder for the download
rawdataDir = 'raw_data'

states = ['LD']
statesRawFolder = 'raw_data'

statesBoundFiles = ['/home/olivierp/jde_coffee/data/Vietnam/LD/aoi/AOI_LD.shp']

# DOWNLOAD parameters
dload = False
# Product to download
product = 'MOD13Q1.006'
# Username for the earthdata website
user = "olivierpfrancois"
# Password for the earthdata website
pwd = "Michele1950"
# Tiles to download.
tiles = ['h28v07']  # ['h28v07']
# Start date for the product download (format YYYY-MM-DD)
#    If None, will default to date of most recent MODIS file on disk if any, 
#    or stop the process
startDownload = None
# End date for the product download (format YYYY-MM-DD)
#    If None, defaults to today
endDownload = None

# MOSAIC parameters
mosaic = False
startMosaic = '2005-01-01'
endMosaic = None

# QUALITY parameters
checkQuality = True

# Get the percentage of missing values for each pixel across history
checkPercMissing = False

maskLand = False

#########################################################################

if dload:
    newHdf = md.downloadMODIS(dstFolder=os.path.join(dst, rawdataDir),
                           pwd=pwd, user=user, tiles=tiles, product=product,
                           startDownload=startDownload, endDownload=endDownload)
else:
    newHdf = []

if mosaic:
    
    thereHdf = [f for f in os.listdir(os.path.join(dst, rawdataDir)) 
                if f.endswith('.hdf')]
    
    if not newHdf and not thereHdf:
        pass
    
    else:
        if not startMosaic:
            print('No start date provided for mosaic. Will mosaic all files on disk')
        print('Starting the mosaic process')
        md.mosaicMODISWrapper(root=dst, srcFolder=os.path.join(dst, rawdataDir),
                       tmpFolder=tempDir, regions=states, regionsOut=statesRawFolder,
                       regionsBoundaries=statesBoundFiles, tiles=tiles,
                       subset='1 0 1 0 0 0 0 0 0 0 0 0', suffix=['NDVI', 'Quality'],
                       startMosaic=startMosaic, endMosaic=endMosaic)

if checkQuality:
    
    # Get all the images on disk
    allNDVI = [os.path.join(dst, 'LD/raw_data', f) for 
               f in os.listdir(os.path.join(dst, 'LD/raw_data')) 
               if f.endswith('NDVI.tif') and '_2005' not in f]
    allNDVI.sort()
    allQuality = [os.path.join(dst, 'LD/raw_data', f) for 
                  f in os.listdir(os.path.join(dst, 'LD/raw_data')) 
                  if f.endswith('Quality.tif') and '_2005' not in f]
    allQuality.sort()
    
    allOut = [f.replace('MODIS/collection6/Vietnam/LD/raw_data', 'Temp') 
              for f in allNDVI]
    
    # Define the dataset
    dataset = zip(allNDVI, allQuality, allOut)
    
    for d in dataset:
        md.maskQualityVI(d[0], d[1], d[2])
    
    '''
    def functionUnpack(args):
        return maskQualityVI(*args)
    
    @context.contextmanager
    def poolcontext(*args, **kwargs):
        pool = mp.Pool(*args, **kwargs)
        yield pool
        pool.terminate()

    with poolcontext(processes=3) as pool:
        pool.map(functionUnpack, dataset)
    
    
    maskQualityVI(dst + '/CO/raw_data/MOD13Q1_2005-01-01_h10v7-8-9_250m_16_days_NDVI.tif',
                  dst + '/CO/raw_data/MOD13Q1_2005-01-01_h10v7-8-9_250m_16_days_Quality.tif',
                  tempDir + '/temp.tif')
    '''

if checkPercMissing:
    # Get the names of all the masked rasters on file
    onDisk = [os.path.join(dst, 'CO/masked_missing', f) for 
                   f in os.listdir(os.path.join(dst, 'CO/masked_missing')) 
                   if f.endswith('NDVI.tif')]
    
    # Dates of these files
    datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
    # Transform into date format
    datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
    
    # Transform into days from start of the year and keep only the unique values
    days = {d.strftime('%j') for d in datesAll}
    # Transform back into a list and order by date
    days = [int(d) for d in days]
    days.sort()
    
    # Get the years for the files on disk
    years = {d.strftime('%Y') for d in datesAll}
    years = [int(y) for y in years]
    years.sort()
    
    for d in days:
        # Get the names of all the rasters for this date
        dates = [datetime(y, 1, 1).date() + timedelta(days=d - 1) for y in years]
        files = [f for f, date in zip(onDisk, datesAll) if date in dates]
        
        # Prepare output name
        outFile = re.sub('([0-9]{4}-[0-9]{2}-[0-9]{2})', str(d), files[0])
        outFile = re.sub('masked', 'perc', outFile)
        
        # Extract the percentage of missing
        md.percMissingStack(images=files, outName=outFile, nodata=None)

'''
# Get the names of all the rasters with percentage missing on file
onDisk = [os.path.join(dst, 'CO/perc_missing', f) for 
                   f in os.listdir(os.path.join(dst, 'CO/perc_missing')) 
                   if f.endswith('.tif')]

# Dates of these files
daysAll = [re.search('MOD13Q1_([0-9]+)_', f).group(1) for f in onDisk]

if maskLand:
    # Get the land mask
    landMask = os.path.join(dst, 'CO/masks/landmask.tif')
    landMask = gdal.Open(landMask)
    
    # Get no data value
    noLand = landMask.GetRasterBand(1).GetNoDataValue()
    
    # Transform into array
    landArray = landMask.GetRasterBand(1).ReadAsArray()
    
    # Loop through the images to mask with the land areas
    for f in onDisk:
        # Import the image
        img = gdal.Open(f)
        
        # Get the no data value
        noImg = img.GetRasterBand(1).GetNoDataValue()
        if not noImg:
            noImg = -99
            img.GetRasterBand(1).SetNoDataValue(noImg)
        
        # Transform into array
        imgArray = img.GetRasterBand(1).ReadAsArray()
        
        # Mask using the land data
        imgArray[landArray == noLand] = noImg
        
        out = md.new_raster_from_base(img, os.path.join(tempDir, os.path.basename(f)),
                                      'GTiff', -99, gdal.GDT_Float32)
        out.GetRasterBand(1).WriteArray(imgArray)
        
        out.FlushCache()
        out = None
        
        # Write values back into file 
        # img.GetRasterBand(1).WriteArray(imgArray)
        
        # img.FlushCache()
        img = None
        
    landMask = None

# Loop through the images and tabulate them
for f, d in zip(onDisk, daysAll):
    outTxt = os.path.join(tempDir, 'miss_' + str(d) + '.txt')
    shp = '/home/olivierp/jde_coffee/data/Colombia/aoi/AOI_cafe_municipios.shp'
    
    img = gdal.Open(f)
    imgArray = img.GetRasterBand(1).ReadAsArray()
    imgArray = imgArray * 1000000
    imgArray.astype('int')
    out = md.new_raster_from_base(img, os.path.join(tempDir, 'temp.tif'),
                                      'GTiff', img.GetRasterBand(1).GetNoDataValue(),
                                      gdal.GDT_UInt32)
    out.GetRasterBand(1).WriteArray(imgArray)
        
    out.FlushCache()
    out = None
    img = None    
    
    img_proc.tabulate_area(regionsFileName=shp,
                           rasterFileName=os.path.join(tempDir, 'temp.tif'),
                           bands=None,
                           polIdFieldName='ID_2', numPix=True,
                           outTxt=outTxt, prefix='val_',
                           numOutDecs=0, alltouch=True)
    
    os.remove(os.path.join(tempDir, 'temp.tif'))
'''
