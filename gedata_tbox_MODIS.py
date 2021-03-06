# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Functions to handle the treatment of MODIS images, 
    Including the download, mosaic, temporal smoothing and averaging.
 """

# Imports

import pymodis as pm
import re, os
#import multiprocessing
import pathos.multiprocessing as mp
from datetime import datetime, timedelta
from osgeo import gdal, ogr
import numpy as np
import functools, math
import gedata_tbox_raster as rt
import mpl_toolkits
mpl_toolkits.__path__.append('/usr/lib/python2.7/dist-packages/mpl_toolkits/')
import matplotlib.pyplot as plt
from matplotlib import cm, dates
from mpl_toolkits.basemap import Basemap, pyproj
import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.patheffects as pe
# import warnings
# warnings.filterwarnings('error')

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# FUNCTIONS


def downloadMODIS(dstFolder, pwd, user, tiles, product, startDownload=None,
                  endDownload=None, satellite='terra'):
    '''
    Function to download modis images from server. Returns a list of the names 
        of the newly downloaded files.
    
    dstFolder (str): Address of the destination folder where to save the 
        downloaded images
    pwd (str): Password for the e4ftl01.cr.usgs.gov website
    user (str): Username for the e4ftl01.cr.usgs.gov website
    tiles (list of str): List of  tiles to download (e.g. 'h13v10')
    product (str): Product to download from server (e.g. 'MOD13Q1.006')
    startDownload (str): Start date for the product download (format YYYY-MM-DD). 
        If None, will default to date of most recent MODIS file on disk if any, 
        or stop the process
    endDownload (str): End date for the product download (format YYYY-MM-DD). 
        If None, defaults to today
    satellite (str): One of 'terra' or 'aqua'
    '''
    
    if not satellite in ['aqua', 'terra']:
        print('MODIS satellite should be aqua or terra')
        return False
    
    if not endDownload:
        # Defaults to today's date
        endDownload = datetime.now().date()
    else:
        endDownload = datetime.strptime(endDownload, '%Y-%m-%d').date()
        
    if startDownload:
        startDownload = datetime.strptime(startDownload, '%Y-%m-%d').date()
    
    # Get the names of the hdf files already downloaded and processed
    # thereTif = [f for f in os.listdir(dstFolder) if f.endswith('.tif')]
    thereHdf = [f for f in os.listdir(dstFolder) if f.endswith('.hdf')]
    # Dates of these files
    if thereHdf:
        datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in thereHdf]
        datesHdf = [datetime.strptime('0101' + d[0:4], "%d%m%Y").date() + timedelta(days=int(d[4:])) 
                     for d in datesHdf]
        # datesTif = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})_', f).group(1) for f in thereTif]
        # datesTif = [datetime.strptime(d, "%Y-%m-%d").date() for d in datesTif]
    else:
        datesHdf = []
        # datesTif = []
    
    if not datesHdf and not startDownload:
        print('There are no files on disk. Choose a start date for the download.')
        return
    
    # Get latest date of the file on disk
    if datesHdf and not startDownload:
        startDownload = max(datesHdf)
        startDownload = min(startDownload, endDownload - timedelta(days=120))
        
    for f, d in zip(thereHdf, datesHdf):
        if d >= startDownload and d <= endDownload:
            try:
                os.remove(os.path.join(dstFolder, f))
                os.remove(os.path.join(dstFolder, f + '.xml'))
            except:
                pass
        
    startDownload = startDownload.strftime('%Y-%m-%d')
        
    endDownload = endDownload.strftime('%Y-%m-%d')
    
    if satellite == 'terra':
        pathData = 'MOLT'
    else:
        pathData = 'MOLA'
    
    # Download 
    down = pm.downmodis.downModis(destinationFolder=dstFolder, password=pwd,
                                  user=user, url="https://e4ftl01.cr.usgs.gov",
                                  tiles=tiles, path=pathData, product=product,
                                  today=startDownload, enddate=endDownload,
                                  jpg=False, debug=False, timeout=30,
                                  checkgdal=True)
    down._connectHTTP()
    down.downloadsAllDay()
    
    # Get the list of downloaded hdf files on disk
    newHdf = [f for f in os.listdir(dstFolder) if f.endswith('.hdf')]
    datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in newHdf]
    datesHdf = [datetime.strptime('0101' + d[0:4], "%d%m%Y").date() + timedelta(days=int(d[4:])) 
                    for d in datesHdf]
    
    newHdf = [f for f, d in zip(newHdf, datesHdf) if 
              d > datetime.strptime(startDownload, '%Y-%m-%d').date()]
    
    if newHdf:
        print(str(len(newHdf)) + ' images downloaded')
        
        # Remove the .txt files
        xml = [f for f in dstFolder if f.endswith('.txt') or f.endswith('.log')]
        for f in xml:
            os.remove(os.path.join(dstFolder, f))
    else:
        print('No new images downloaded')
    
    return newHdf

   
def mosaicMODISWrapper(root, srcFolder, tmpFolder, regions, regionsOut,
                       regionsBoundaries, tiles, subset, suffix, nodataOut=None,
                       startMosaic=None, endMosaic=None):
    '''
    Function to mosaic the tiles of modis images and clip them to a series of 
        regions. Does not return anything.
    
    root (str): Address of root folder where  regions folders are located
    srcFolder (str): Full address where the input downloaded tiles are located
    tmpFolder (str): Full address of folder where the temporary mosaic before 
        cutting to regions will be stored. 
        !!!!! THIS FOLDER SHOULD BE EMPTY TO START, IT WILL BE EMPTIED AT THE END
    regions (list of str): Names of the regions to process. Each region should 
        have a folder inside root with the same name.  
    regionsOut (str): Name of the folder inside the regions folders where the 
        mosaiced and clipped images should be stored. 
        It should be the same for all the regions
    regionsBoundaries (list of str): List the full address of the shapefiles 
        with the regions boundaries
    tiles (list of str): List of  tiles to mosaic (e.g. 'h13v10')
    subset (str): string with the layers to extract from the hdf images 
        (e.g. '1 0 0 0 0 0 0 0 0 0 0 0')
    suffix (list): list of string, suffixes to assign for naming each of the 
        subsets
    nodataOut (list): None or list of the same length as suffix, with the no 
        data value to use for each layer extracted 
        when masking by the regions. Can leave one no data value as None to 
        use the no data value of the raster
    startMosaic (str): Starting date for the files to mosaic. If None, will 
        process all the files found on the disk.
    endMosaic (str): Ending date for the files to mosaic. If None, 
        defaults to today
    '''
    
    if not endMosaic:
        endMosaic = datetime.now()
        endMosaic = endMosaic.date()
    else:
        endMosaic = datetime.strptime(endMosaic, '%Y-%m-%d').date()
    
    if startMosaic:
        startMosaic = datetime.strptime(startMosaic, '%Y-%m-%d').date()
        
    # Get the files on disk
    thereHdf = [f for f in os.listdir(srcFolder) if f.endswith('.hdf')]
    # Dates of these files
    datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in thereHdf]
    datesHdfD = [datetime.strptime('0101' + d[0:4], "%d%m%Y").date() + timedelta(days=int(d[4:]) - 1) 
                 for d in datesHdf]
    
    # Get the dates to process
    if startMosaic:
        dates = {d for d, D in zip(datesHdf, datesHdfD) 
                 if D >= startMosaic and D <= endMosaic}
    else:
        dates = {d for d, D in zip(datesHdf, datesHdfD) if D <= endMosaic}
    dates = sorted(list(dates))
    
    # Get the tiles for the output name
    tilesH = [re.sub('v[0-9]+', '', t) for t in tiles]
    tilesH = [re.findall(r'\d+', t) for t in tilesH]
    tilesH = {int(t[0]) for t in tilesH}
    tilesH = list(tilesH)
    tilesH.sort()
    tilesV = [re.sub('h[0-9]+', '', t) for t in tiles]
    tilesV = [re.findall(r'\d+', t) for t in tilesV]
    tilesV = {int(t[0]) for t in tilesV}
    tilesV = list(tilesV)
    tilesV.sort()
    
    # Work by date
    for d in dates:
        # Transform into actual date
        dlong = datetime.strptime(d, '%Y%j').date()
        
        # Progress
        print('Processing date ' + dlong.strftime('%Y-%m-%d'))
        
        # Get the file names for that date
        files = [f for f in thereHdf if 'A' + d in f]
        
        # Check if all the tiles were completed
        complete = True
        for t in tiles:
            if not complete:
                break
            
            complete = any(t in f for f in files)
        
        if not complete:
            print('Not all tiles were downloaded for date ' + 
                  dlong.strftime('%Y-%m-%d') + ". Cannot process.")
            continue
        
        # Create the output name for the mosaic file
        outName = re.sub('\..+', '', files[0])  # Get the product name
        outName = (outName + '_' + dlong.strftime('%Y-%m-%d') + '_h' + 
                   '-'.join(map(str, tilesH)) + 'v' + 
                   '-'.join(map(str, tilesV)) + '_250m_16_days')
        
        mosaicMODIS(images=[os.path.join(srcFolder, f) for f in files],
                    subset=subset,
                    suffixes=suffix,
                    tempFolder=tmpFolder,
                    outFile='/'.join([tmpFolder, outName + '.tif']))
        
        if not nodataOut:
            nodataOut = [None for _ in suffix]
            
        # Loop through the mosaics
        for l in range(len(suffix)):
            
            # Loop through the regions to clip and mask the mosaics
            for s, shp in zip(regions, regionsBoundaries):
                
                # Clip the mosaic by the extent of the the shapefile
                clipMaskRasterByShp(shp=shp,
                                    raster='/'.join([tmpFolder, outName + '_' + 
                                                     suffix[l] + '.tif']),
                                    outRaster=root + '/' + s + '/' + 
                                    regionsOut + '/' + outName + '_' + 
                                    suffix[l] + '.tif',
                                    clipR=True,
                                    maskR=False,
                                    nodataOut=nodataOut[l])
                
            # Remove intermediary file
            os.remove('/'.join([tmpFolder, outName + '_' + suffix[l] + '.tif']))
        
        # Remove the .xml, .hdf and .vrt intermediary files
        xml = [f for f in os.listdir(os.path.join(tmpFolder)) if f.endswith('.xml') or 
                    f.endswith('.hdf') or f.endswith('.vrt')]
        for f in xml:
            os.remove(os.path.join(tmpFolder, f))


def mosaicMODIS(images, subset, suffixes, tempFolder, outFile):
    '''
    Creates mosaic from a set of MODIS images. It creates one composite 
        image per layer selected in the subset.
    
    images (list): list of addresses of images to mosaic. They should be from 
        the same dates but from different tiles.
    subset (str): string with the layers to extract from the hdf images 
        (e.g. '1 0 0 0 0 0 0 0 0 0 0 0')
    suffix (list): list of string, suffixes to assign for naming each of the subsets
    tmpFolder (str): Full address of folder where the temporary mosaic before 
        cutting to regions will be stored. 
        !!!!! THIS FOLDER SHOULD BE EMPTY TO START, IT WILL BE EMPTIED AT THE END
    outFile (str): Full name of the output file. Should be complete with 
        extension .tif. The final output will have the additional suffix for 
        each subset
    '''
    
    # Export each subset layer separately
    layers = [i for i, x in enumerate(list(subset)) if x == '1']
    
    for i in range(len(layers)):
        # Make a copy of the lists
        layersCopy = list(layers)
        subsetCopy = list(subset)
        # Remove the layer being considered
        layersCopy.pop(i)
        if type(layersCopy) is int:
            layersCopy = [layersCopy]
        # Replace all the others by 0 and return into a string
        for l in list(layersCopy):
            subsetCopy[l] = '0'
        subsetCopy = ''.join(subsetCopy)
        
        # Mosaic the files and transform to tif
        mos = pm.convertmodis_gdal.createMosaicGDAL(hdfnames=images,
                subset=subsetCopy, outformat='GTiff')  # outformat='HDF4Image'
        mos.run(output=os.path.join(tempFolder, 'temp_' + suffixes[i] + '.tif'))  # '.hdf'
        # mos.write_mosaic_xml(os.path.join(tmpFolder,'outName'))
        # mos.write_vrt(os.path.join(tmpFolder,outName))
    
        # Change the projection to lat/long
        inRaster = '/'.join([tempFolder, 'temp_' + suffixes[i] + '.tif'])
        outRaster = outFile.replace('.tif', '_' + suffixes[i] + '.tif')
        
        cmd = 'gdalwarp -overwrite -t_srs EPSG:4326 -r near -of GTiff -co BIGTIFF=IF_NEEDED %s %s' % (inRaster, outRaster)
        os.system(cmd)
        
        os.remove(inRaster)


def clipMaskRasterByShp(shp, raster, outRaster, clipR=True, maskR=True,
                        dataToMask=None, nodataOut=None, alltouch=False):
    '''
    Function clips a raster by the extent of a shapefile and masks 
        the areas outside of the polygons as an option
    
    shp (str): Full address of the shapefile to clip
    raster (str): Full address of the raster to be clipped.
    outRaster (str): Full address of the output raster
    clipR (bool): Whether the image should be clipped by the shapefile
    maskR (bool): Whether the image should be masked by the shapefile
    dataToMask (list): Optional list of values of the data to mask
        outside of the polygon. 
        If None, all the values will be masked
    nodataOut (num): Optional no data value for the output raster. If None, 
        will use the no data value from the input raster.
    alltouch (bool): True/False if all the pixels touched by the polygons 
        should be counted as covered by the polygons.
    '''
    
    if not clipR and not maskR:
        return False
    
    # Open the data source and read in the extent
    driver = ogr.GetDriverByName('ESRI Shapefile')
    maskDS = driver.Open(shp, 0)  # 0 means read-only. 1 means writeable.
    maskLayer = maskDS.GetLayer()
    minX, maxX, minY, maxY = maskLayer.GetExtent()
    
    # read the image
    image = gdal.Open(raster)
    
    # Get the geotransform and the projection
    geoTrans = image.GetGeoTransform()
    projection = image.GetProjection()
    # Get number of columns and rows
    mosaicCols = image.RasterXSize
    mosaicRows = image.RasterYSize
    
    # Get the pixel position of the shapefile corners
    ulX, ulY = rt.world2Pixel(geoTrans, minX, maxY)
    lrX, lrY = rt.world2Pixel(geoTrans, maxX, minY)
    
    # Create a new geoTransform for the output raster
    geoRegion = list(geoTrans)
    if ulX >= 0:
        # Update the origin if the image is larger than the shapefile
        geoRegion[0] = minX
    if ulY >= 0:
        # Update the origin if the image is larger than the shapefile
        geoRegion[3] = maxY
    
    # Check the pixel position to correct if limits are outside the 
    # image
    ulX = min(max(ulX, 0), mosaicCols - 1)
    lrX = min(max(lrX, 0), mosaicCols - 1)
    ulY = min(max(ulY, 0), mosaicRows - 1)
    lrY = min(max(lrY, 0), mosaicRows - 1)
    
    if ulX == lrX or ulY == lrY:
        print('Shapefile and raster don\'t overlap')
        return False
    
    # Get the band
    band = image.GetRasterBand(1)
    
    # Get the data type
    dType = band.DataType
    
    nodataImg = nodataOut
    if not nodataImg:
        nodataImg = band.GetNoDataValue()
    
    if not nodataImg:
        print('No NODATA information in raster and none provided')
        return False
    
    # Extract the image information in that extent if clipping or get all
    if clipR:
        # using readAsArray xoffset, yoffset, xextent, yextent
        clip = band.ReadAsArray(ulX, ulY, int(lrX - ulX), int(lrY - ulY)).astype(np.float)
    else:
        clip = band.ReadAsArray().astype(np.float)

    # Remove rasters to save
    if os.path.isfile(outRaster):
        os.remove(outRaster)
    
    # Create an empty raster with the right size
    rows, cols = clip.shape
    driver = gdal.GetDriverByName("GTiff")
    new_raster = driver.Create(outRaster, cols, rows, 1, getattr(gdal, rt.getGDALTypeFromNumber(dType)),
                               options=['COMPRESS=LZW', 'BIGTIFF=IF_NEEDED'])
    new_raster.SetProjection(projection)
    new_raster.SetGeoTransform(geoRegion)

    new_raster.GetRasterBand(1).SetNoDataValue(nodataImg)
    new_raster.GetRasterBand(1).Fill(nodataImg)
    
    if maskR:
        # Rasterize the region into that raster
        if alltouch:
            gdal.RasterizeLayer(new_raster, [1], maskLayer, burn_values=[1.],
                                options=['ALL_TOUCHED=TRUE'])
        else:
            gdal.RasterizeLayer(new_raster, [1], maskLayer, burn_values=[1.],
                                options=['ALL_TOUCHED=FALSE'])
        
        # Get the resulting raster band as an array
        new = new_raster.GetRasterBand(1).ReadAsArray().astype(np.float)
        
        # Mask the values outside the shapefile, either totally or selectively
        if dataToMask:
            for d in dataToMask:
                clip[(new == nodataImg) & (clip == d)] = nodataImg
        else:
            clip[new == nodataImg] = nodataImg
    
    image = None
    
    # Export the values
    new_raster.GetRasterBand(1).WriteArray(clip, 0, 0)
    
    # Close the raster
    new_raster.FlushCache()
    new_raster = None
    
    # Crop the mosaic and mask the raster
    # inRaster = tmpFolder + '/' + outName
    # outRaster = root + '/' + s + '/' + regionsOut + '/' + outName
    # Call gdalwarp using the os system command
    # cmd = 'gdalwarp -q -cutline %s -crop_to_cutline %s %s' % (shp, inRaster, outRaster)
    # os.system(cmd)
    # Call gdalwarp using the subprocess command
    # out = subprocess.call(['/usr/bin/gdalwarp', '-q', '-overwrite', '-cutline', shp, '-crop_to_cutline', 
    #                       '-of GTIFF', '-co BIGTIFF=IF_NEEDED', inRaster, outRaster])
    # print(out)


def smoothMODISWrapper(root, regions, regionsIn, regionsOut, startSmooth, endSmooth, regWindow, avgWindow,
                startSaveSmooth=None, endSaveSmooth=None, algo='Swets', parallel=False, nCores=None):
    '''
    Function to do a temporal smoothing of a time series of identical images. 
        Does not return anything, saves smoothed images on disk.
    The input files should contain the date in their name in the format 
        '_([0-9]{4}-[0-9]{2}-[0-9]{2})'
    The processed files will have the same name plus a 'smooth_' prefix
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region should 
        have a folder inside root with the same name.
    regionsIn (str): Name of the folder inside the regions folders where to 
        find the input files to be smoothed 
            It should be the same for all the regions 
    regionsOut (str): Name of the folder inside the regions folders where to 
        save the smoothed images 
        It should be the same for all the regions
    startSmooth (str): Starting date for the files to mosaic. If None, will 
        process all the files found on the disk.
    endSmooth (str): Ending date for the files to mosaic. If None, defaults 
        to today
    regWindow (int): size of the regression window (see Swets et al. 
        for details)
    avgWindow (int): sie of the averaging window (see Swets et al. for details)
    startSaveSmooth (str): Starting date for the files to save. If None, will 
        save all the processed files
    endSaveSmooth (str): Ending date for the files to save. If None, will save 
        all the processed files
    algo(str): Name of algorithm to use for the smoothing. Can be Swets or 
        Savitzky
    '''
    
    if not endSmooth:
        endSmooth = datetime.now()
        endSmooth = endSmooth.date()
    else:
        endSmooth = datetime.strptime(endSmooth, '%Y-%m-%d').date()
    if not startSmooth:
        startSmooth = endSmooth - timedelta(days=365)
    else:
        startSmooth = datetime.strptime(startSmooth, '%Y-%m-%d').date()
    
    # Transform into date format
    if not startSaveSmooth:
        startSaveSmooth = endSmooth - timedelta(days=183)
    else:
        startSaveSmooth = datetime.strptime(startSaveSmooth, '%Y-%m-%d').date()
    if not endSaveSmooth:
        endSaveSmooth = endSmooth
    else:
        endSaveSmooth = datetime.strptime(endSaveSmooth, '%Y-%m-%d').date()
    
    # Declare parallel workers if option
    if parallel and len(regions) > 1:
        if not nCores:
            nCores = mp.cpu_count()
        
        p = mp.Pool(nCores)
    
        pp = functools.partial(smoothRegionWrap, rootR=root,
                               regionsIn=regionsIn,
                               regionsOut=regionsOut,
                               startSmooth=startSmooth,
                               endSmooth=endSmooth,
                               regWindow=regWindow,
                               avgWindow=avgWindow,
                               startSaveSmooth=startSaveSmooth,
                               endSaveSmooth=endSaveSmooth,
                               algo=algo)
                
        p.map(pp, regions)
        
        # Close the threads
        p.close()
        p.join()
        
    else:
        # Loop through the regions to do the smoothing for each
        for r in regions:
            smoothRegionWrap(rootR=root, r=r,
                             regionsIn=regionsIn,
                             regionsOut=regionsOut,
                             startSmooth=startSmooth,
                             endSmooth=endSmooth,
                             regWindow=regWindow,
                             avgWindow=avgWindow,
                             startSaveSmooth=startSaveSmooth,
                             endSaveSmooth=endSaveSmooth,
                             algo=algo)


def smoothRegionWrap (r, rootR, regionsIn, regionsOut,
                      startSmooth, endSmooth, regWindow, avgWindow,
                      startSaveSmooth, endSaveSmooth, algo):
    
    print('Processing region ' + str(r) + '...')
        
    # Import all the raw modis images on disk
    onDisk = [os.path.join(rootR, r, regionsIn, f) 
              for f in os.listdir(os.path.join(rootR, r, regionsIn)) 
              if f.endswith('.tif') and 'NDVI' in f]
    
    # Dates of these files
    datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                for f in onDisk]
    # Transform into date format
    datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
    
    # Keep only the files and dates within the dates to process
    onDisk = [f for f, d in zip(onDisk, datesAll) 
              if d >= startSmooth and d <= endSmooth]
    datesAll = [d for d in datesAll if d >= startSmooth and d <= endSmooth]
    
    # Sort the two list by date
    datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
    
    # Create a mask to know which of the dates to save
    toSave = [d >= startSaveSmooth and d <= endSaveSmooth for d in datesAll]
    
    # Smooth the images
    smoothSeries(inRasters=onDisk, toSave=toSave,
                 outFolder=rootR + '/' + r + '/' + regionsOut,
                 regWindow=regWindow,
                 avgWindow=avgWindow,
                 algo=algo, blockXSize=256, blockYSize=256)


def smoothSeries(inRasters, toSave, outFolder, regWindow, avgWindow,
                 algo='Swets', blockXSize=256, blockYSize=256):
    '''
    Implements the Swets or Savitzky temporal smoothing over a series of images.
    The images should be provided in the correct temporal order
    
    The processed files will have the same name plus a 'smooth_' prefix
    
    inRasters (list): Full addresses of the rasters to be smoothed.
    toSave (list): List of booleans of the same size as inRasters. 
        For each raster, whether the smoothed version should be saved to disk 
        or not.
    outFolder (str): Full address of the folder where to export the smoothed 
        rasters.
    regWindow (int): size of the regression window (see Swets et al. 
        for details)
    avgWindow (int): sie of the averaging window (see Swets et al. for details)
    algo(str): Name of algorithm to use for the smoothing. Can be Swets or 
        Savitzky
    blockXSize (int): X size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    blockYSize (int): Y size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    '''
    
    if not len(inRasters) == len(toSave):
        print('toSave should be a list of boolean with the same length as the inputs')
        return False
    
    # Import all the images to process in the right order
    toProcess = [gdal.Open(f) for f in inRasters]
    
    # Get the no data value if any
    nodata = toProcess[0].GetRasterBand(1).GetNoDataValue()
    if not nodata:
        nodata = np.nan
    
    # Remove rasters to save if any
    for s, f in zip(toSave, inRasters):
        if s and os.path.isfile(outFolder + '/smooth_' + os.path.basename(f)):
            os.remove(outFolder + '/smooth_' + os.path.basename(f))
    
    # Create empty copies of these files to use for the smoothed data only for the files to save
    processed = [rt.newRasterFromBase(p, outFolder + '/smooth_' + os.path.basename(f),
                                      'GTiff', nodata, gdal.GDT_Float32) 
                    if s else False for s, p, f in zip(toSave, toProcess, inRasters)]
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = toProcess[1].GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Get the number of blocks in x and y directions based on block size
    xBlocks = int(round(xsize / blockXSize)) + 1
    yBlocks = int(round(ysize / blockYSize)) + 1
    
    totBlocks = xBlocks * yBlocks
    progress = 1
    
    for xStep in range(xBlocks):
        for yStep in range(yBlocks):
            print('     Processing block ' + str(progress) + ' of ' + 
                  str(totBlocks))
            progress += 1
            
            blocks = [rt.readRasterBlock(p, xStep * blockXSize, yStep * blockYSize,
                                      blockXSize, blockYSize) 
                                      for p in toProcess]
            
            # Bring the blocks together into one single array
            blocks = np.dstack(blocks)
            # Recast the type
            blocks = blocks.astype(np.float32)
            
            # Do the smoothing for each pixel
            if algo == 'Swets':
                blocks = smoothingSwets(blocks, regWindow, avgWindow, nodata)
            elif algo == 'Savitzky':
                blocks = smoothingSavitzky(blocks, nodata)
            
            # Return to regular ndvi values between -1 and 1
            blocks[blocks != nodata] = np.divide(blocks[blocks != nodata], 10000.)
            
            # Change the values in the output raster
            blocks = np.dsplit(blocks, len(toProcess))
            for s, p, b in zip(toSave, processed, blocks):
                if s:
                    p.GetRasterBand(1).WriteArray(b[:, :, 0], xStep * blockXSize,
                                                  yStep * blockYSize)
    
    # Close the rasters
    for s, p in zip(toSave, processed):
        if s:
            p.FlushCache()
            p = None
    
    for p in zip(toSave, toProcess):            
        p = None


def smoothingSwets(block, regWindow, avgWindow, nodata=None, minLim=None, maxLim=None):
    '''
    Function to do the temporal smoothing of a 3D numpy array.
    The function will loop through all the columns in the first two dimensions, 
        smoothing the data in the 3rd dimension.
    The method is that of Swets et al. (1999, A weighted least-squares approach 
        to temporal NDVI smoothing. ASPRS Annual Conference).
    It returns an array of the same shape with smoothed values, 
        where the no data value is np.nan if not provided.
        The function assumes that numpy nan is the no data value if it is not 
        provided.
    
    The Swets et al. method runs a moving (weighted) regression window along 
        the values in the input pixels. It then
    averages the predicted values using a second moving window the get the 
        final smoothed values.
    
    block (numpy array): Array to be smoothed
    regWindow (int): size of the regression window (see Swets et al. 
        for details)
    avgWindow (int): sie of the averaging window (see Swets et al. for details)
    nodata (int): value for the no data in the input array if any. If None, 
        assumed to be np.nan. The output array will be returned with np.nan as 
        the no data value if nodata is not provided.
    minLim (num): minimum value that the pixel can take. 
        The output will be limited to that value
    maxLim (num): maximum value that the pixel can take. 
        The output will be limited to that value
    '''
    
    extent = block.shape
    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Check if entire pixel vector is nan, just return it
            if (nodata and np.all(pixel == nodata)):
                pixel.fill(nodata)
                block[X, Y, :] = pixel
                continue
            elif np.all(np.isnan(pixel)):
                pixel.fill(np.nan)
                block[X, Y, :] = pixel
                continue
            
            # Get the shape of the original data
            originShape = pixel.shape
            
            # Change the data type
            pixel = pixel.astype(np.float32)
            
            # Replace no data values by nan if needed
            if nodata:
                pixel[pixel == nodata] = np.nan
            
            # Keep track of the nodata value positions
            miss = np.isnan(pixel)
            
            # Reshape the data
            pixel = np.reshape(pixel, (len(pixel), 1))
            
            # Interpolate missing values if any just for the
            # sake of the smoothing. The no data will be 
            # added back to the output.
            # Get a mask of the nan values
            if np.isnan(pixel).any():
                pixel = np.reshape(pixel, len(pixel))
                nans = np.isnan(pixel)
                # Get the actual indices for these values and the non nan values
                notNans = ~nans
                notNans = notNans.nonzero()[0]
                nans = nans.nonzero()[0]
                # Interpolate
                pixel[nans] = np.interp(nans, notNans, pixel[notNans])
                # Reshape
                pixel = np.reshape(pixel, (len(pixel), 1))
            
            # Create a vector of weights for the regression
            weights = np.copy(pixel)
            weights[0, 0] = 0.5  # Assume first value is a middle point
            for i in range(1, len(pixel) - 1, 1):
                # Local high values get a weight of 1.5, local middle values 
                # get a weight of 0.5 and local low values only 0.005
                if (pixel[i - 1, 0] < pixel[i, 0] and 
                    pixel[i, 0] > pixel[i + 1, 0]):
                    weights[i, 0] = 1.5
                elif ((pixel[i - 1, 0] <= pixel[i, 0] and 
                       pixel[i, 0] <= pixel[i + 1, 0]) or 
                       (pixel[i - 1, 0] >= pixel[i, 0] and 
                        pixel[i, 0] >= pixel[i + 1, 0])):
                    weights[i, 0] = 0.5
                elif (pixel[i - 1, 0] > pixel[i, 0] and 
                      pixel[i, 0] < pixel[i + 1, 0]):
                    weights[i, 0] = 0.005
            # For the last point
            if (pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 2, 0] and 
                pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 3, 0]):
                # If the last data point is greater than the previous 2, then 
                # assume it is a high point
                weights[len(pixel) - 1, 0] = 1.5
            elif (pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 2, 0] and 
                  pixel[len(pixel) - 1, 0] < pixel[len(pixel) - 3, 0]):
                # If the last data point is greater than the previous only, 
                # then assume it is a middle point
                weights[len(pixel) - 1, 0] = 0.5
            elif (pixel[len(pixel) - 1, 0] < pixel[len(pixel) - 2, 0] and 
                  pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 3, 0] and 
                  weights[len(pixel) - 3, 0] != 0.005):
                # If less than the previous but more than the one before that 
                # is not a low point
                weights[len(pixel) - 1, 0] = 0.5
            else:
                # If the last point is less than the last 2, or greater than 
                # only one of the two, assume it is a low point
                weights[len(pixel) - 1, 0] = 0.002
            
            # Create a matrix with the data for this pixel and for the weights
            # For the data:
            # Each column will be the same
            dataRaw = np.repeat(pixel, len(pixel) + regWindow, axis=1)
            # Set to nan all data for each column besides the data used for the 
            # regressions
            # Lower triangle indices below the diagonal
            ltri = np.tril_indices(n=len(pixel), k=-1, m=len(pixel) + regWindow)
            dataRaw[ltri] = np.nan
            utri = np.triu_indices(n=len(pixel), k=0, m=len(pixel))
            dataRaw[:, (regWindow):][utri] = np.nan
            # Remove the first two and last 3 columns, since they don't have 
            # enough points for a regression
            dataRaw = dataRaw[:, 2:(len(pixel) + regWindow - 3)]
            
            # For the weights:
            weights = np.repeat(weights, len(pixel) + regWindow - 5, axis=1)
            weights[np.isnan(dataRaw)] = np.nan
            
            # Create an empty array for the results of the regressions
            dataReg = np.zeros(dataRaw.shape)
            
            # Estimate the regression for each column of dataRaw
            for i in range(len(pixel) + regWindow - 5):
                # Prepare regression data
                y = dataRaw[:, i][~np.isnan(dataRaw[:, i])]  # dependent
                w = weights[:, i][~np.isnan(dataRaw[:, i])]  # weights
                x = range(1, len(y) + 1)  # independent
                x = np.asarray(x)
                x = x.astype(np.float32)
                x = np.reshape(x, y.shape)
                # Estimate potential outliers
                # Will be an outlier if greater than 3 (only extreme outliers 
                # should be picked up)
                out = np.divide(np.absolute(y - np.mean(y)), np.std(y))
                out = np.reshape(out, y.shape)
                # Remove outliers before regression
                yout = y[out < 3]
                wout = w[out < 3]
                xout = x[out < 3]
                # Compute parameters of regression
                numerator = (np.sum(wout) * 
                             np.sum(np.multiply(np.multiply(wout, xout), yout)) - 
                             np.sum(np.multiply(wout, yout)) * 
                             np.sum(np.multiply(wout, xout)))
                denominator = (np.sum(wout) * 
                               np.sum(np.multiply(wout, np.square(xout))) - 
                               np.square(np.sum(np.multiply(wout, xout))))
                b = np.divide(numerator, denominator)
                numerator = (np.sum(np.multiply(wout, yout)) - b * 
                             np.sum(np.multiply(wout, xout)))
                denominator = np.sum(wout)
                a = np.divide(numerator, denominator)
                # Compute the predicted values from the regression
                dataReg[:, i][~np.isnan(dataRaw[:, i])] = a + b * x
            
            dataReg[np.isnan(dataRaw)] = np.nan
            
            # Combination of the results from the regression for each point
            # Now we combine for each point the results from the moving regression window
            # We take into account the results from avg.window regression windows, centered around the point, unless we are at the edges
            smoothed = []  # Will hold the smoothed results
            t = int(np.floor(avgWindow / 2.))  # number of predicted values from regressions to take into account around the center
            for i in range(len(pixel)):
                x = dataReg[i, :][~np.isnan(dataReg[i, :])]
                if i < np.floor(regWindow / 2.) and len(x) < regWindow:
                    center = int(np.floor(regWindow / 2.)) - (regWindow - len(x))
                    res = np.mean(x[max(0, center - t):max(1, center + t + 1)])
                elif i == len(pixel) - 1:
                    res = np.mean(
                        x[(int(np.floor(regWindow / 2.)) - t - 1):(int(np.floor(regWindow / 2.)) + t + 1)])
                else:
                    res = np.mean(
                        x[(int(np.floor(regWindow / 2.)) - t):(int(np.floor(regWindow / 2.)) + t + 1)])
                smoothed.append(res)
            
            smoothed = np.asarray(smoothed)
            smoothed = np.reshape(smoothed, originShape)
            
            if minLim:
                smoothed[~np.isnan(smoothed) & smoothed < minLim] = minLim
            if maxLim:
                smoothed[~np.isnan(smoothed) & smoothed > maxLim] = maxLim
            
            if nodata:
                smoothed[np.isnan(smoothed)] = nodata
            
            # Replace the original missing pixels back
            if nodata:
                smoothed[miss] = nodata
            else:
                smoothed[miss] = np.nan
            
            block[X, Y, :] = smoothed
    
    return block


def smoothingSavitzky(block, nodata=None):
    '''
    Function to do the temporal smoothing of a 3D numpy array.
    The function will loop through all the columns in the first two dimensions, 
        smoothing the data in the 3rd dimension.
    The method is that of Savitzky-Golay, using a window of 7 and a 2nd order 
        polynomial.
    It returns an array of the same shape with smoothed values, 
        where the no data value is np.nan if not provided.
        The function assumes that numpy nan is the no data value if it is not 
        provided.
    
    The Savitzky-Golay method runs a moving (weighted) window along the 
        values in the input pixels. It simply 
    averages the values in the window to predict the output value.
    A window of 7 is still used for edages, but different weights are applied.
    
    block (numpy array): Array to be smoothed
    nodata (int): value for the no data in the input array if any. If None, 
        assumed to be np.nan. The output array will be returned with np.nan as 
        the no data value if nodata is not provided.
    '''
    
    extent = block.shape
    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Check if entire pixel vector is nan, just return it
            if (nodata and np.all(pixel == nodata)):
                pixel.fill(nodata)
                block[X, Y, :] = pixel
                continue
            elif np.all(np.isnan(pixel)):
                pixel.fill(np.nan)
                block[X, Y, :] = pixel
                continue
            
            # Get the shape of the original data
            originShape = pixel.shape
            
            # Reshape the data
            pixel = np.reshape(pixel, len(pixel))
            pixel = pixel.astype(np.float32)
            
            # Replace no data values by nan if needed
            if nodata:
                pixel[pixel == nodata] = np.nan
            
            # Get a copy of the pixel with only the non nan values
            toSmooth = pixel[~np.isnan(pixel)]
            
            # The moving window of 7 doesn't work if there aren't enough data
            if len(toSmooth) < 7:
                pixel.fill(np.nan)
                block[X, Y, :] = pixel
                continue
            
            for _ in range(3):
                # Take a first smoothing of the pixel
                smoothed = savitzkyAvg(toSmooth)
                # Replace the values that are lower in the original pixel by 
                # the values in the original smoothed pixel
                np.copyto(toSmooth, smoothed, where=toSmooth < smoothed)
            
            # Replace the smoothed values in the original pixel
            pixel[~np.isnan(pixel)] = toSmooth
            
            # Interpolate missing values if any
            # Get a mask of the nan values
            if np.isnan(pixel).any():
                nans = np.isnan(pixel)
                # Get the actual indices for these values and the non nan values
                notNans = ~nans
                notNans = notNans.nonzero()[0]
                nans = nans.nonzero()[0]
                # Interpolate
                pixel[nans] = np.interp(nans, notNans, pixel[notNans])
                
            # Reshape
            pixel = np.reshape(pixel, originShape)
            
            if nodata:
                smoothed[np.isnan(smoothed)] = nodata
            
            # Replace the smoothed pixel in the block
            block[X, Y, :] = pixel
    
    return block


def savitzkyAvg(pixel):
    '''
    This function averages a multi-dates pixel using a window 7 for the 
        Savitzky Golay and returns a pixel of the same shape
    
    pixel (numpy array): Numpy vector to be averaged along using Savitzky Golay.
    '''
    
    out = np.copy(pixel)
    
    for i in range(len(pixel)):
        if i < 3:
            # Take the moving window
            avg = pixel[0:7]
            # Do the weighted sum
            if i == 0:
                out[i] = np.sum(
                    np.multiply(avg,
                                np.array([0.1190, -0.0714, -0.1429, -0.0952, 0.0714, 0.3571, 0.7619])))
            elif i == 1:
                out[i] = np.sum(
                    np.multiply(avg,
                                np.array([-0.0714, -0.0000, 0.0714, 0.1429, 0.2143, 0.2857, 0.3571])))
            elif i == 2:
                out[i] = np.sum(
                    np.multiply(avg,
                                np.array([-0.1429, 0.0714, 0.2143, 0.2857, 0.2857, 0.2143, 0.0714])))
        elif i > (len(pixel) - 4):
            # Take the moving window
            avg = pixel[(len(pixel) - 7):]
            # Do the weighted sum
            if i == len(pixel) - 3:
                out[i] = np.sum(
                    np.multiply(avg,
                                np.array([0.0714, 0.2143, 0.2857, 0.2857, 0.2143, 0.0714, -0.1429])))
            elif i == len(pixel) - 2:
                out[i] = np.sum(
                    np.multiply(avg,
                                np.array([0.3571, 0.2857, 0.2143, 0.1429, 0.0714, -0.0000, -0.0714])))
            elif i == len(pixel) - 1:
                out[i] = np.sum(
                    np.multiply(avg,
                                np.array([0.7619, 0.3571, 0.0714, -0.0952, -0.1429, -0.0714, 0.1190])))
        else:
            # Take the moving window
            avg = pixel[i - 3:i + 4]
            # Do the weighted sum
            out[i] = np.sum(np.multiply(avg, np.array([-0.0952, 0.1429, 0.2857, 0.3333, 0.2857, 0.1429, -0.0952])))                
        
    return(out)
            

def createBaseline(root, regions, varieties, regionsIn, regionsOut, startRef,
                   endRef, mask=None, outModelRaster=None, parallel=False, nCores=None):
    '''
    Function to create baseline images with the decile values over the period 
        of interest. 
    Does not return anything, saves baselines images on disk (10 layer rasters).
    The input files should come from the smoothing process of modis data.
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region 
        should have a folder inside root with the same name.
    varieties (list of lists): varieties to consider for each region. 
        The function will create a series of baselines for each regions and 
        each varieties for that region. The name of the variety will be used 
        to name the baselines.
        Elements in the list of lists are string, for example 
        [['arabica','robusta'],['arabica'],...].
    regionsIn (str): Name of the folder inside the regions folders where to 
        find the input modis rasters to be used for the baselines. 
        These files are outputs of the smoothing. 
        It should be the same for all the regions.
    regionsOut (str): Name of the folder inside the regions folders where to 
        save the baselines images.
        It should be the same for all the regions.
    startRef (str): Starting year for the files to consider for the baseline 
        (included).
    endRef (str): Ending year for the files to consider for the baseline 
        (included).
    mask (list of lists): None or list of lists of full addresses of the masks 
        (.tif) to use for each of the regions and varieties.
        Each list element of the overall list should have the same length as 
        the corresponding list for varieties or be None/empty.
        The masks should be at the same resolution and stack perfectly with 
        the modis images.
        If None, no masking will be done.
    outModelRaster (list of lists): None or list of lists of full addresses of 
        the rasters (.tif) to use as model for the output.
        or each of the regions and varieties.
    '''
    
    if parallel:
        if not nCores:
            nCores = mp.cpu_count()
        
        p = mp.Pool(nCores)
    
    # Loop through the regions
    for r in range(len(regions)):
        # Get the names of all the smoothed rasters on file
        onDisk = [os.path.join(root, regions[r], regionsIn, f) 
                  for f in os.listdir(os.path.join(root, regions[r], regionsIn)) 
                  if f.endswith('.tif')]
        
        # Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                    for f in onDisk]
        # Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        
        # Transform into days from start of the year and keep only the unique values
        days = {d.strftime('%j') for d in datesAll}
        # Transform back into a list and order by date
        days = [int(d) for d in days]
        days.sort()
        
        for v in range(len(varieties[r])):
            
            # Check if there is a mask and a model for the output raster
            if mask and mask[r] and mask[r][v]:
                maskD = mask[r][v]
            else:
                maskD = None
            if outModelRaster and outModelRaster[r] and outModelRaster[r][v]:
                outModelRasterD = outModelRaster[r][v]
            else:
                outModelRasterD = None
            
            # Loop through the dates to create a baseline raster for each
            if parallel:
                pp = functools.partial(decileWrap,
                                       startRef=startRef,
                                       endRef=endRef,
                                       onDisk=onDisk,
                                       datesAll=datesAll,
                                       root=root,
                                       region=regions[r],
                                       regionsOut=regionsOut,
                                       variety=varieties[r][v],
                                       maskD=maskD,
                                       outModelRasterD=outModelRasterD)
                        
                p.map(pp, days)
                
            else:
                for d in days:
                    decileWrap(d=d, startRef=startRef,
                               endRef=endRef,
                               onDisk=onDisk,
                               datesAll=datesAll,
                               root=root,
                               region=regions[r],
                               regionsOut=regionsOut,
                               variety=varieties[r][v],
                               maskD=maskD,
                               outModelRasterD=outModelRasterD)
        
    if parallel:
        # Close the threads
        p.close()
        p.join()


def decileWrap(d, startRef, endRef, onDisk, datesAll,
               root, region, regionsOut, variety, maskD,
               outModelRasterD):
    # Get the names of all the rasters for this date
    dates = [datetime(y, 1, 1).date() + timedelta(days=d - 1) for 
             y in range(startRef, endRef + 1, 1)]
    # dates = [datetime.strptime('0101'+str(y), '%d%m%Y').date()+timedelta(days=d-1) for y in range(startRef,endRef+1,1)]
    files = [f for f, date in zip(onDisk, datesAll) if date in dates]
    
    # Prepare outName
    outName = (root + '/' + region + '/' + 
               regionsOut + '/' + 'ndvi_deciles_0to100pct_ref_period-' + 
               str(startRef) + '-' + str(endRef) + 
               '_day' + str(d) + '_date-' + 
               dates[0].strftime('%b-%d') + '_' + 
               variety + '.tif')
    
    # Create decile raster
    createDecileRaster(images=files,
                  outFile=outName,
                  mask=maskD,
                  outModelRaster=outModelRasterD,
                  blockXSize=256, blockYSize=256)

                
def createDecileRaster(images, outFile, mask=None, outModelRaster=None, blockXSize=256, blockYSize=256):
    '''
    Takes a set of images that perfectly overlap (from different dates) and 
    create an output image with 10 layers, each of which has the value of 
    a decile for each pixel.
    
    images (list):Full addresses of the images to be used as reference for the 
        calculation of the deciles.
    outFile (str): Full address of the output decile raster.
    mask (str): None or full address of the mask to be used over the images for 
        where to calculate the decile references.
    outModelRaster (str): None or full address of the raster to use 
        as model for the output. In that case the output raster will 
        be warped using that raster as model before export.
    blockXSize (int): X size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    blockYSize (int): Y size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    '''
    
    # Import all the images to use for estimating the deciles
    sourceImg = [gdal.Open(f) for f in images]
    
    # Get the no data value
    nodata = sourceImg[0].GetRasterBand(1).GetNoDataValue()
    
    # Create an empty copy in memory
    toProcess = [rt.newRasterFromBase(p, '', 'MEM',
                                     nodata, gdal.GDT_Float32, bands=1) for
                 p, f in zip(sourceImg, images)]
    
    # Fill it with the values
    for p, s in zip(toProcess, sourceImg):
        p.GetRasterBand(1).WriteArray(s.GetRasterBand(1).ReadAsArray())
        
        s = None
    
    # Import mask if there is one and apply it
    if mask:
        # Import the mask raster as an array
        n = gdal.Open(mask)
        nArray = n.GetRasterBand(1).ReadAsArray()
        #nodataMask = n.GetRasterBand(1).GetNoDataValue()
        
        # Create the mask
        nArray = np.logical_or(
            np.logical_or(
                nArray <= 0, np.isnan(nArray)),
                                nArray > 1)
        
        '''
        nArray = np.logical_or(
            np.logical_or(
                nArray == 0, np.isnan(nArray)),
                                nArray == nodataMask)
        '''
        
        # Apply the mask
        for p in toProcess:
            pArray = p.GetRasterBand(1).ReadAsArray()
            pArray[nArray] = nodata
            p.GetRasterBand(1).WriteArray(pArray)
        
        n = None
        nArray = None
        pArray = None
        
    # Change resolution if there is a model
    if outModelRaster:
        '''
        #Create an empty raster in memory with the wanted resolution
        base = gdal.Open(os.path.join(root,regions[r],regionsIn,files[0]))
        projection = base.GetProjection()
        geotransform = base.GetGeoTransform()
        extentX = base.rasterXSize*geotransform[1]
        extentY = base.regions[r]asterYSIze*geotransform[5]
        #Change the resolution to what is wanted
        geotransform[1] = outResolution[0]
        geotransform[5] = outResolution[1]
        ncols = int(round(extentX/outResolution[0]))
        nrows = int(round(extentY/outResolution[1]))
        #Create the empty raster
        driver = gdal.GetDriverByName('MEM')
        templateRes = driver.Create('', ncols, nrows, 1, gdal.GDT_Byte)
        templateRes.SetGeoTransform(geotransform)
        templateRes.SetProjection(projection)
        base = None
        '''
        
        # Change the resolution of the rasters
        outModel = gdal.Open(outModelRaster)
        toProcess = [rt.warpRaster(p, outModel, resampleOption='average', outputURI=None, outFormat='MEM') 
                     for p in toProcess]
        
    # Remove existing raster if any
    if os.path.isfile(outFile):
        os.remove(outFile)
    
    # Create an empty copy with 10 layers to use for storing the deciles
    processed = rt.newRasterFromBase(toProcess[0], outFile, 'GTiff',
                                     nodata, gdal.GDT_Float32, bands=10)
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = toProcess[1].GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Get the number of blocks in x and y directions based on block size
    xBlocks = int(round(xsize / blockXSize)) + 1
    yBlocks = int(round(ysize / blockYSize)) + 1
    
    for xStep in range(xBlocks):
        for yStep in range(yBlocks):
            
            block = [rt.readRasterBlock(p, xStep * blockXSize, yStep * blockYSize, blockXSize, blockYSize) for p in toProcess]
            
            # Bring the blocks together into one single array
            block = np.dstack(block)
            # Recast the type
            block = block.astype(np.float32)
            
            # Estimate the deciles for each pixel
            deciles = estimateDeciles(block, nodata)
            
            # Change the values in the output raster
            deciles = np.dsplit(deciles, 10)
            for i in range(10):
                processed.GetRasterBand(i + 1).WriteArray(deciles[i][:, :, 0], xStep * blockXSize, yStep * blockYSize)
                    
    # Close the rasters
    for p in toProcess:
        p.FlushCache()
        p = None
    processed.FlushCache()
    processed = None


def estimateDeciles(block, nodata):
    '''
    Takes a block of data (numpy array) with 3 dimensions and estimates for
    each pixel in the third dimension value of the deciles.
    the function returns a block with the same number of pixels and 10 layers, 
    one for each decile.
    
    block (array): 3D numpy array with the pixels to process
    nodata (num): value to be considered as no data when processing the pixels.
    '''
    
    extent = block.shape
    deciles = np.empty((extent[0], extent[1], 10))

    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Replace the no data value by nan
            pixel[pixel == nodata] = np.nan
            
            # Get the number of nan values in the pixel
            nbNodata = np.sum(np.isnan(pixel))
            
            # Return no data if there are less than 6 valid values
            if len(pixel) - nbNodata < 6:
                deciles[X, Y, :].fill(nodata)
                continue
            
            # Remove the nan values and reshape the data
            pixel = pixel[~np.isnan(pixel)]
            pixel = np.reshape(pixel, (len(pixel), 1))
            pixel = pixel.astype(np.float32)
            
            # Find the deciles
            deciles[X, Y, :] = np.reshape(np.percentile(pixel, np.arange(0, 100, 10)), deciles[X, Y, :].shape)
    
    return deciles


def rankDatesDeciles(root, regions, varieties, regionsIn, refDecilesIn, startRank, endRank, mask, minDensity=None):
    '''
    Function to rank each pixel of an image based on the baseline for the 
        same day of year.
    Does not return anything, saves ranked image on disk.
    The input files should come from the smoothing and baseline functions.
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region should 
        have a folder inside root with the same name.
    varieties (list of lists): varieties to consider for each region. 
        The function will create a series of baselines for 
        each regions and each varieties for that region. The name of the 
        variety will be used to name the baselines.
        Elements in the list of lists are string, for example 
        [['arabica','robusta'],['arabica'],...]
    regionsIn (str): Name of the folder inside the regions folders where to 
        find the input modis rasters to be used for the ranking. These files 
        are outputs of the smoothing. It should be the same for all the regions 
    refDecilesIn (str): Name of the folder inside the regions folders where 
        to find the input baseline rasters to be used for the ranking. 
        These files are outputs of the baseline creation. 
        It should be the same for all the regions 
    regionsOut (str): Name of the folder inside the regions folders where to 
        save the baselines images 
        It should be the same for all the regions
    startRank (str): Starting date for the files to rank (included).
    endRank (str): Ending date for the files to rank (included).
    mask (list of lists): None or list of lists of full paths of the masks 
        (.tif) to use for each of the regions and varieties.
        Each list element of the overall list should have the same length as 
        the corresponding list for varieties or be None/empty.
        The masks should be at the same resolution and stack perfectly with 
        the modis images.
        The mask should contain for each pixel either nodata or the density of 
        the crop masked in the pixel (0 < d < 1).
        If None, no masking will be done.
    minDensity (decimal): None or list of values between 0 and 1. Minimum 
        density of crop of interest to filter out pixels in the ranking. 
        All pixels with less than these densities will be masked, with one 
        output per density.
        None by default (no masking).
    '''
    
    # Transform into date format
    if not endRank:
        endRank = datetime.now()
        endRank = endRank.date()
    else:
        endRank = datetime.strptime(endRank, '%Y-%m-%d').date()
    if not startRank:
        startRank = endRank - timedelta(days=60)
    else:
        startRank = datetime.strptime(startRank, '%Y-%m-%d').date()
    
    # Loop through the regions to do the ranking for each
    for r in range(len(regions)):
        print('Ranking region ' + str(regions[r]) + '...')
        
        # Import all the smooth modis images on disk
        onDisk = [f for f in os.listdir(os.path.join(root, regions[r], regionsIn)) if f.endswith('.tif')]
        
        # Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        # Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
        # Keep only the files and dates within the dates to process
        onDisk = [f for f, d in zip(onDisk, datesAll) if d >= startRank and d <= endRank]
        datesAll = [d for d in datesAll if d >= startRank and d <= endRank]
        # Sort the two list by date
        # datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
        # Transform into days from start of the year and keep only the unique values
        days = [int(d.strftime('%j')) for d in datesAll]
        
        # Get the names of the decile baselines on disk
        baseFiles = [f for f in os.listdir(os.path.join(root, regions[r], refDecilesIn)) if f.endswith('.tif')]
        
        for v in range(len(varieties[r])):
            
            # Loop through the days to rank each date
            for day, date, fileS in zip(days, datesAll, onDisk):
                
                # Check if there is a mask and a model for the output raster
                if mask and mask[r] and mask[r][v]:
                    maskD = mask[r][v]
                else:
                    maskD = None
                
                # Get the address of the decile raster
                baseFile = [f for f in baseFiles if 'day' + str(day) + '_' in f and varieties[r][v] in f]
                baseFile = os.path.join(root, regions[r], refDecilesIn, baseFile[0])
                
                # Create the output name
                outName = os.path.join(root + '/' + regions[r] + '/' + 
                                         'ndvi_' + date.strftime('%Y-%m-%d') + 
                                         '_CompareToDecile_0BelowMin_110AboveMax_' + 
                                         varieties[r][v] + '.tif')
                
                # Rank the pixels of the image
                estimateRankRaster(image=os.path.join(root, regions[r], regionsIn, fileS),
                                   deciles=baseFile,
                                   densityMask=maskD,
                                   outFile=outName,
                                   minDensity=minDensity,
                                   blockXSize=256, blockYSize=256)


def estimateRankRaster(image, deciles, outFile, densityMask=None,
                       minDensity=None, blockXSize=256, blockYSize=256):
    '''
    Function takes an image and a raster with the reference deciles for that 
    image (from estimateDeciles) and returns for each pixel its position
    compared to the decile values using estimateRank. 
    The image and the reference deciles raster can have different resolutions, 
    in which case teh resolution of th eimage will be matched to the resolution 
    of the deciles raster.
    
    image (str): Full address of the image to be ranked.
    deciles (str): Full address of the image containing the reference 
        deciles for that image.
    densityMask (str): None or mask to be used for where to compare the image 
        to the deciles based on the density information in the mask. 
        Only the pixels with a minimum density of minDensity will be ranked.
    outFile (str): Full address of the output file.
    minDensity (list): None or list of density values to be used as threshold 
        for ranking. There will be one output per minimum value
    blockXSize (int): X size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    blockYSize (int): Y size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    '''
    
    if densityMask:
        # Import the mask
        nanMask = gdal.Open(densityMask)
    
    # Import the image to rank
    smoothImg = gdal.Open(image)
    
    # Get no data value
    nodata = smoothImg.GetRasterBand(1).GetNoDataValue()
    
    # Import the decile baseline
    baseImg = gdal.Open(deciles)
    
    # Adapt the resolution of the smooth images to the baseline if needed
    geoSmooth = smoothImg.GetGeoTransform()
    geoBase = baseImg.GetGeoTransform()
    reproject = [1 for a, b in zip(geoSmooth, geoBase) if not a == b]
    if reproject:
        smoothImg = rt.warpRaster(smoothImg, baseImg, resampleOption='average', outputURI=None, outFormat='MEM')
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = baseImg.GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Get the number of blocks in x and y directions based on the block size
    xBlocks = int(round(xsize / blockXSize)) + 1
    yBlocks = int(round(ysize / blockYSize)) + 1
    band = None
    
    if not densityMask:
        minDensity = [0]
    
    else:
        if not minDensity:
            minDensity = [0]
    
    for m in minDensity:
        if m >= 1 or m < 0:
            continue
        
        # Prepare output name
        outRaster = re.sub('.tif', '_maskedbelow' + str(int(m * 100)) + 'pct.tif', outFile)
        
        # Remove existing raster if any
        if os.path.isfile(outRaster):
            os.remove(outRaster)
        
        # Create an empty copy for storing the deciles comparisons
        processed = rt.newRasterFromBase(baseImg, outRaster, 'GTiff',
                                         nodata, gdal.GDT_Int16, bands=1)
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                # Read the block from the images
                blockSmooth = rt.readRasterBlock(smoothImg, xStep * blockXSize,
                                              yStep * blockYSize,
                                              blockXSize, blockYSize)
                
                blockBase = [rt.readRasterBlock(baseImg, xStep * blockXSize,
                                             yStep * blockYSize,
                                             blockXSize, blockYSize, band=b + 1) 
                             for b in range(baseImg.RasterCount)]
                
                if densityMask:
                    # Read the block from the mask 
                    blockMask = rt.readRasterBlock(nanMask, xStep * blockXSize,
                                                yStep * blockYSize,
                                                blockXSize, blockYSize)
                
                    if m > 0:
                        blockMask = blockMask < m
                    elif m == 0:
                        blockMask = blockMask <= m
                    
                    # Apply the mask
                    blockSmooth[blockMask] = nodata
                    for b in blockBase:
                        b[blockMask] = nodata
                    
                # Bring the blocks together into one single array
                blockBase = np.dstack(blockBase)
                
                # Recast the type to be sure
                blockSmooth = blockSmooth.astype(np.float32)
                blockBase = blockBase.astype(np.float32)
            
                # Estimate the placement for each pixel
                ranks = estimateRank(block=blockSmooth, ref=blockBase, nodata=nodata)
                
                # Change the values in the output raster
                processed.GetRasterBand(1).WriteArray(ranks, xStep * blockXSize, yStep * blockYSize)
    
    # Close the rasters
    smoothImg = None
    baseImg = None
    processed.FlushCache()
    processed = None
    
    if densityMask:
        nanMask = None


def estimateRank(block, ref, nodata):
    '''
    Takes a block of data (array) from an image and rank each pixel against 
    a block of reference data with the deciles. 
    The pixel is ranked as -3, -2, -3, 10, 20, ..., 100, and 110. 
    110 is if the pixel is above the max decile value
    10 the pixel is at the min value exactly, -1 it is between the min and 
    0.75 of the mean, -2 between 0.75 and 0.5 of the mean, and finally -3 if
    less tha 0.5 the min.
    
    block (array): array of pixels to be ranked
    ref (array): array of reference decile values, should have same dimensions 
        as block but with 10 layers in third dimension
    nodata (num): value to be considered as no data, so not ranked.
    '''
    
    extent = block.shape
    ranking = np.empty((extent[0], extent[1]), dtype=int)
    
    for X in range(extent[0]):
        for Y in range(extent[1]):
            testPixel = np.copy(block[X, Y])
            refPixel = np.copy(ref[X, Y, :])
            refPixel = np.reshape(refPixel, (len(refPixel), 1))
            
            # Replace the no data value by nan
            refPixel[np.logical_or(refPixel == nodata, np.isinf(refPixel))] = np.nan
            testPixel[np.logical_or(testPixel == nodata, np.isinf(testPixel))] = np.nan
            
            # return no data if the test pixel is nan or the ref pixel is nan
            if np.isnan(testPixel) or np.all(np.isnan(refPixel)):
                ranking[X, Y] = nodata
                continue
            
            # Get the position of the pixel in the reference deciles
            refPixel = np.reshape(refPixel, len(refPixel))
            
            rank = np.searchsorted(refPixel, testPixel) * 10
            
            # If the value is above the deciles, returns length of the vector
            if len(refPixel) * 10 == rank:
                rank = 110
            # If the value is at min or below, returns 0
            if rank == 0:
                if testPixel == refPixel[0]:
                    rank = 10
                elif testPixel / refPixel[0] >= 0.75:
                    rank = -1
                elif testPixel / refPixel[0] >= 0.5:
                    rank = -2
                else:
                    rank = -3
            
            ranking[X, Y] = rank
    
    return np.reshape(ranking, extent)


def avgRegionRasterWrap(regionIn, avgWeights, weightField=None,
                        startAvg=None, endAvg=None, alltouch=False):
    '''
    Function to compute the average ndvi value for a region of interest.
    It returns the dictionary with the averages.
    This is a wrapper for avgRegionRaster, simply pulling out the 
    images between the dates requested.
    
    regionIn (str): Full address of the folder inside the regions folders where to 
        find the input modis rasters to be used for the averaging. These 
        files are outputs of the smoothing. 
    avgWeights (str): Raster (single layer) or shapefile with the weights 
        to use for the averages. If shapefile, will be rasterized. 
        If raster, will be resampled to the correct resolution if needed. 
        In each case the information in each pixel/polygon should be a density 
        of crop of interest (0 < d < 1).
    weightField (str): Only used if avgWeights is a shapefile. Name of the 
        field containing the weight/density information.
    startAvg (str): Starting date for the files to consider in the averaging 
        (included). There will be one value per date.
    endAvg (str): Ending date for the files to consider in the averaging 
        (included).
    alltouch (boolean): true or false, whether all the pixels touched by the 
        region should be considered in the average or only the pixels with 
        their centroid inside the region.
        This parameter is used when rasterizing the shapefile with the 
        weights and will therefore only be used if avgWeights is a shapefile.
    '''
    
    # Transform into date format
    if not startAvg:
        startAvg = datetime.now() - timedelta(days=365)
        startAvg = startAvg.date()
    else:
        startAvg = datetime.strptime(startAvg, '%Y-%m-%d').date()
    if not endAvg:
        endAvg = datetime.now()
        endAvg = endAvg.date()
    else:
        endAvg = datetime.strptime(endAvg, '%Y-%m-%d').date()
        
    # Import all the smooth modis images on disk
    onDisk = [os.path.join(regionIn, f) 
              for f in os.listdir(regionIn) if f.endswith('.tif')]
    # Dates of these files
    datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                for f in onDisk]
    # Transform into date format
    datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
    # Keep only the files and dates within the dates to process
    onDisk = [f for f, d in zip(onDisk, datesAll) if d >= startAvg and d <= endAvg]
    datesAll = [d for d in datesAll if d >= startAvg and d <= endAvg]
    
    if not onDisk:
        return {}
    
    # Compute the averages for all the dates
    avgRegion = avgRegionRaster(images=onDisk,
                                datesImg=datesAll,
                                weightsRaster=avgWeights,
                                weightField=weightField,
                                alltouch=alltouch,
                                blockXSize=256, blockYSize=256)
    
    if not avgRegion:
        return False
    
    return avgRegion


def avgRegionRaster(images, datesImg, weightsRaster=None, weightField=None,
                    alltouch=False, blockXSize=256, blockYSize=256):
    '''
    Takes a set of images and computes the average using the weights in 
    weightsRaster if any for each of them.
    Returns a dictionary with the average value for each image. The dictionary 
    keys are the dates of the images as provided in datesImg
    
    images (list): list of full addresses of images to average.
    datesImg (list): list of the same length as images, providing for each
        the date of the image (or some unique index). They should be unique.
    weightsRaster (str): Raster (single layer) or shapefile with the weights 
        to use for the averages. If shapefile, will be rasterized. 
        If raster, will be resampled to the correct resolution if needed. 
        In each case the information in each pixel/polygon should be a density 
        of crop of interest (0 < d < 1).
    weightField (str): Only used if avgWeights is a shapefile. Name of the 
        field containing the weight/density information.
    alltouch (boolean): true or false, whether all the pixels touched by the 
        region should be considered in the average or only the pixels with 
        their centroid inside the region.
        This parameter is used when rasterizing the shapefile with the 
        weights and will therefore only be used if avgWeights is a shapefile.
    blockXSize (int): X size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    blockYSize (int): Y size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    '''
    
    if weightsRaster:
        # Get a base image as a reference for format
        baseImg = gdal.Open(images[0])
        
        # Rasterize the gridded weights if needed or simply import it
        if weightsRaster.endswith('.shp'):
            if not weightField:
                print('The name of the field with the densities needs to be specified to average')
                return(False)
            
            # Import the vector layer
            # Open the shapefile
            driver = ogr.GetDriverByName('ESRI Shapefile')
            dataSource = driver.Open(weightsRaster, 0)  # 0 means read-only. 1 means writeable.
            
            # Create layer
            inVecLayer = dataSource.GetLayer(0)
            
            # Prepare an empty raster to rasterize the shapefile
            weightsRaster = rt.newRasterFromBase(baseImg, 'temp', 'MEM', -1, gdal.GDT_Float32) 
            
            # Transform alltouch
            if alltouch:
                alltouch = 'TRUE'
            else:
                alltouch = 'FALSE'
            
            # Rasterize the vector layer:
            gdal.RasterizeLayer(weightsRaster, [1], inVecLayer,
                                options=['ALL_TOUCHED=' + alltouch, 'ATTRIBUTE=' + weightField])
            
            inVecLayer = None
            
        elif weightsRaster.endswith(('.tif', '.TIF')):
            weightsRaster = gdal.Open(weightsRaster)
            
            # Change the resolution of the raster to match the images if needed
            
            geoSmooth = weightsRaster.GetGeoTransform()
            geoBase = baseImg.GetGeoTransform()
            reproject = [1 for a, b in zip(geoSmooth, geoBase) if not a == b]
            if reproject:
                weightsRaster = rt.warpRaster(weightsRaster, baseImg, resampleOption='nearest', outputURI=None, outFormat='MEM')
        
        baseImg = None
        nodataWeights = weightsRaster.GetRasterBand(1).GetNoDataValue()
    
    # Prepare an empty dictionary to hold the results
    average = {}
    
    # Loop through the images to compute the average for each
    for img, date in zip(images, datesImg):
        
        # Import the image
        baseImg = gdal.Open(img)
        
        # Loop through the blocks to compute the average for each
        
        # Get the size of the rasters to identify the limits of the blocks to loop through
        band = baseImg.GetRasterBand(1)
        # Get the no data value
        nodataBase = band.GetNoDataValue()
        # Get the size of the raster
        xsize = band.XSize
        ysize = band.YSize
        # Get the number of blocks in x and y directions based on the block size
        xBlocks = int(round(xsize / blockXSize)) + 1
        yBlocks = int(round(ysize / blockYSize)) + 1
        band = None
        
        sumNdvi = []
        sumWeights = []
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                # Read the block from the image
                blockBase = rt.readRasterBlock(baseImg,
                                            xStep * blockXSize, yStep * blockYSize,
                                            blockXSize, blockYSize)
                
                # Recast the type to be sure
                blockBase = blockBase.astype(np.float32)
                
                # Replace the no data values by 0
                blockBase[np.logical_or(
                    np.logical_or(blockBase == nodataBase, np.isnan(blockBase)),
                    np.logical_or(blockBase > 1, blockBase < -1))] = 0.
                
                if weightsRaster:
                    # Read the block from the weights
                    blockWeight = rt.readRasterBlock(weightsRaster,
                                                  xStep * blockXSize, yStep * blockYSize,
                                                  blockXSize, blockYSize)
                    
                    # Recast the type to be sure
                    blockWeight = blockWeight.astype(np.float32)
                
                    # Replace the no data values by 0
                    blockWeight[np.logical_or(
                        np.logical_or(blockWeight == nodataWeights, np.isnan(blockWeight)),
                        np.logical_or(blockWeight > 1, blockWeight < 0))] = 0.
                    
                    # Set the no data values in the image as zero weight
                    blockWeight[blockBase == 0.] = 0.
                    
                else:
                    blockWeight = np.copy(blockBase)
                    
                    # Replace the non 0 data values by 1
                    blockWeight[blockWeight != 0.] = 1.
                
                # Estimate the weighted sum for each pixel
                sumNdvi.append(np.sum(np.multiply(blockBase, blockWeight)))
                sumWeights.append(np.sum(blockWeight))
                
        # Combine for the entire image
        sumNdvi = np.sum(sumNdvi) / np.sum(sumWeights)
        
        # Add to the output dictionary along with the date
        average[date.strftime('%Y-%m-%d')] = sumNdvi
    
        # Close the raster
        baseImg = None
    
    return(average)


def computeQualityIndexNdviWrap(regionIn, avgWeights, weightField=None,
                                missingValue=None, startAvg=None, endAvg=None,
                                alltouch=False):
    # Transform into date format
    if not startAvg:
        startAvg = datetime.now() - timedelta(days=365)
        startAvg = startAvg.date()
    else:
        startAvg = datetime.strptime(startAvg, '%Y-%m-%d').date()
    if not endAvg:
        endAvg = datetime.now()
        endAvg = endAvg.date()
    else:
        endAvg = datetime.strptime(endAvg, '%Y-%m-%d').date()
        
    # Import all the modis images on disk
    onDisk = [os.path.join(regionIn, f) 
              for f in os.listdir(regionIn) if f.endswith('.tif')]
    # Dates of these files
    datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                for f in onDisk]
    # Transform into date format
    datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
    # Keep only the files and dates within the dates to process
    onDisk = [f for f, d in zip(onDisk, datesAll) if d >= startAvg and d <= endAvg]
    datesAll = [d for d in datesAll if d >= startAvg and d <= endAvg]
    
    if not onDisk:
        return {}
    
    # Compute the averages for all the dates
    avgRegion = computeQualityIndexNdvi(images=onDisk,
                                        datesImg=datesAll,
                                        missingValue=missingValue,
                                        weightsRaster=avgWeights,
                                        weightField=weightField,
                                        alltouch=alltouch,
                                        blockXSize=256, blockYSize=256)
    
    if not avgRegion:
        return False
    
    return avgRegion
    

def computeQualityIndexNdvi(images, datesImg, missingValue=None,
                            weightsRaster=None, weightField=None,
                    alltouch=False, blockXSize=256, blockYSize=256):
    
    if weightsRaster:
        # Get a base image as a reference for format
        baseImg = gdal.Open(images[0])
        
        # Rasterize the gridded weights if needed or simply import it
        if weightsRaster.endswith('.shp'):
            if not weightField:
                print('The name of the field with the densities needs to be specified to average')
                return(False)
            
            # Import the vector layer
            # Open the shapefile
            driver = ogr.GetDriverByName('ESRI Shapefile')
            dataSource = driver.Open(weightsRaster, 0)  # 0 means read-only. 1 means writeable.
            
            # Create layer
            inVecLayer = dataSource.GetLayer(0)
            
            # Prepare an empty raster to rasterize the shapefile
            weightsRaster = rt.newRasterFromBase(baseImg, 'temp', 'MEM', -1, gdal.GDT_Float32) 
            
            # Transform alltouch
            if alltouch:
                alltouch = 'TRUE'
            else:
                alltouch = 'FALSE'
            
            # Rasterize the vector layer:
            gdal.RasterizeLayer(weightsRaster, [1], inVecLayer,
                                options=['ALL_TOUCHED=' + alltouch, 'ATTRIBUTE=' + weightField])
            
            inVecLayer = None
            
        elif weightsRaster.endswith(('.tif', '.TIF')):
            weightsRaster = gdal.Open(weightsRaster)
            
            # Change the resolution of the raster to match the images if needed
            
            geoSmooth = weightsRaster.GetGeoTransform()
            geoBase = baseImg.GetGeoTransform()
            reproject = [1 for a, b in zip(geoSmooth, geoBase) if not a == b]
            if reproject:
                weightsRaster = rt.warpRaster(weightsRaster, baseImg,
                                            resampleOption='nearest',
                                            outputURI=None, outFormat='MEM')
        
        baseImg = None
        nodataWeights = weightsRaster.GetRasterBand(1).GetNoDataValue()
    
    # Prepare an empty dictionary to hold the results
    quality = {}
    
    # Loop through the images to compute the average for each
    for img, date in zip(images, datesImg):
        
        # Import the image
        baseImg = gdal.Open(img)
        
        # Loop through the blocks to compute the average for each
        
        # Get the size of the rasters to identify the limits of the blocks to loop through
        band = baseImg.GetRasterBand(1)
        # Get the no data value
        nodataBase = band.GetNoDataValue()
        
        # Update the value to consider as missing if not provided
        if not missingValue:
            missingValue = nodataBase
            
        # Get the size of the raster
        xsize = band.XSize
        ysize = band.YSize
        # Get the number of blocks in x and y directions based on the block size
        xBlocks = int(round(xsize / blockXSize)) + 1
        yBlocks = int(round(ysize / blockYSize)) + 1
        band = None
        
        sumWeights = []
        sumAll = []
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                # Read the block from the image
                blockBase = rt.readRasterBlock(baseImg,
                                            xStep * blockXSize, yStep * blockYSize,
                                            blockXSize, blockYSize)
                
                # Recast the type to be sure
                blockBase = blockBase.astype(np.float32)
                
                # Replace the no data values by 0
                # blockBase[np.logical_or(
                #    np.logical_or(blockBase == nodataBase, np.isnan(blockBase)),
                #    np.logical_or(blockBase > 1, blockBase < -1))] = 0.
                
                # Create a mask of the block to only keep 
                # the missing value
                baseMask = blockBase == missingValue
                    
                if weightsRaster:
                    # Read the block from the weights
                    blockWeight = rt.readRasterBlock(weightsRaster,
                                                  xStep * blockXSize, yStep * blockYSize,
                                                  blockXSize, blockYSize)
                    
                    # Recast the type to be sure
                    blockWeight = blockWeight.astype(np.float32)
                
                    # Replace the no data values by 0
                    blockWeight[np.logical_or(
                        np.logical_or(blockWeight == nodataWeights, np.isnan(blockWeight)),
                        np.logical_or(blockWeight > 1, blockWeight < 0))] = 0.
                    
                else:
                    # Create an array of same size as block with only ones
                    blockWeight = np.ones(blockBase.shape)
                    
                # Set the no data values in the image as zero weight
                # blockWeight[blockBase == 0.] = 0.
                
                # Estimate the weighted sum for each pixel
                sumWeights.append(np.sum(np.multiply(baseMask, blockWeight)))
                sumAll.append(np.sum(blockWeight))       
        
        # Combine for the entire image
        sumWeights = np.sum(sumWeights) / np.sum(sumAll)
        
        # Add to the output dictionary along with the date
        quality[date.strftime('%Y-%m-%d')] = sumWeights
    
        # Close the raster
        baseImg = None
    
    return(quality)


def avgRegionQualWrap(regionIn, maskedIn, avgWeights, weightField=None,
                      startAvg=None, endAvg=None, alltouch=False):
    '''
    ALTERNATIVE FUNCTION TO AVGREGIONRASTERWRAP TO CALCULTATE THE 
    AVERAGES ONLY ON GOOD QUALITY PIXELS
    
    Function to compute the average ndvi value for a region of interest.
    It returns the dictionary with the averages.
    This is a wrapper for avgRegionRaster, simply pulling out the 
    images between the dates requested.
    
    regionIn (str): Full address of the folder inside the regions folders where to 
        find the input modis rasters to be used for the averaging. These 
        files are outputs of the smoothing. 
    avgWeights (str): Raster (single layer) or shapefile with the weights 
        to use for the averages. If shapefile, will be rasterized. 
        If raster, will be resampled to the correct resolution if needed. 
        In each case the information in each pixel/polygon should be a density 
        of crop of interest (0 < d < 1).
    weightField (str): Only used if avgWeights is a shapefile. Name of the 
        field containing the weight/density information.
    startAvg (str): Starting date for the files to consider in the averaging 
        (included). There will be one value per date.
    endAvg (str): Ending date for the files to consider in the averaging 
        (included).
    alltouch (boolean): true or false, whether all the pixels touched by the 
        region should be considered in the average or only the pixels with 
        their centroid inside the region.
        This parameter is used when rasterizing the shapefile with the 
        weights and will therefore only be used if avgWeights is a shapefile.
    '''
    
    # Transform into date format
    if not startAvg:
        startAvg = datetime.now() - timedelta(days=365)
        startAvg = startAvg.date()
    else:
        startAvg = datetime.strptime(startAvg, '%Y-%m-%d').date()
    if not endAvg:
        endAvg = datetime.now()
        endAvg = endAvg.date()
    else:
        endAvg = datetime.strptime(endAvg, '%Y-%m-%d').date()
        
    # Import all the smooth modis images on disk
    onDisk = [os.path.join(regionIn, f) 
              for f in os.listdir(regionIn) if f.endswith('.tif')]
    # Dates of these files
    datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                for f in onDisk]
    # Get the corresponding quality files to only take into account 
    # good quality pixels
    qualAll = [os.path.join(maskedIn, f) 
                for f in os.listdir(maskedIn) if f.endswith('.tif')]
    qualDisk = []
    for d in datesAll:
        add = [f for f in qualAll if d in f]
        if not add:
            add = ['']
        qualDisk = qualDisk + add 
    qualAll = None
    # Transform into date format
    datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
    # Keep only the files and dates within the dates to process
    onDisk = [f for f, d in zip(onDisk, datesAll) if d >= startAvg and d <= endAvg]
    datesAll = [d for d in datesAll if d >= startAvg and d <= endAvg]
    
    if not onDisk:
        return {}
    
    # Compute the averages for all the dates
    avgRegion = avgRegionQual(images=onDisk,
                                datesImg=datesAll,
                                maskedQual=qualDisk,
                                weightsRaster=avgWeights,
                                weightField=weightField,
                                alltouch=alltouch,
                                blockXSize=256, blockYSize=256)
    
    if not avgRegion:
        return False
    
    return avgRegion


def avgRegionQual(images, datesImg, maskedQual, weightsRaster=None, weightField=None,
                    alltouch=False, blockXSize=256, blockYSize=256):
    '''
    Takes a set of images and computes the average using the weights in 
    weightsRaster if any for each of them.
    Returns a dictionary with the average value for each image. The dictionary 
    keys are the dates of the images as provided in datesImg
    
    images (list): list of full addresses of images to average.
    datesImg (list): list of the same length as images, providing for each
        the date of the image (or some unique index). They should be unique.
    weightsRaster (str): Raster (single layer) or shapefile with the weights 
        to use for the averages. If shapefile, will be rasterized. 
        If raster, will be resampled to the correct resolution if needed. 
        In each case the information in each pixel/polygon should be a density 
        of crop of interest (0 < d < 1).
    weightField (str): Only used if avgWeights is a shapefile. Name of the 
        field containing the weight/density information.
    alltouch (boolean): true or false, whether all the pixels touched by the 
        region should be considered in the average or only the pixels with 
        their centroid inside the region.
        This parameter is used when rasterizing the shapefile with the 
        weights and will therefore only be used if avgWeights is a shapefile.
    blockXSize (int): X size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    blockYSize (int): Y size of the block to be processed at once from the 
        rasters (code proceeds by block for memory management).
    '''
    
    if weightsRaster:
        # Get a base image as a reference for format
        baseImg = gdal.Open(images[0])
        
        # Rasterize the gridded weights if needed or simply import it
        if weightsRaster.endswith('.shp'):
            if not weightField:
                print('The name of the field with the densities needs to be specified to average')
                return(False)
            
            # Import the vector layer
            # Open the shapefile
            driver = ogr.GetDriverByName('ESRI Shapefile')
            dataSource = driver.Open(weightsRaster, 0)  # 0 means read-only. 1 means writeable.
            
            # Create layer
            inVecLayer = dataSource.GetLayer(0)
            
            # Prepare an empty raster to rasterize the shapefile
            weightsRaster = rt.newRasterFromBase(baseImg, 'temp', 'MEM', -1,
                                                 gdal.GDT_Float32) 
            
            # Transform alltouch
            if alltouch:
                alltouch = 'TRUE'
            else:
                alltouch = 'FALSE'
            
            # Rasterize the vector layer:
            gdal.RasterizeLayer(weightsRaster, [1], inVecLayer,
                                options=['ALL_TOUCHED=' + alltouch,
                                         'ATTRIBUTE=' + weightField])
            
            inVecLayer = None
            
        elif weightsRaster.endswith(('.tif', '.TIF')):
            weightsRaster = gdal.Open(weightsRaster)
            
            # Change the resolution of the raster to match the images if needed
            
            geoSmooth = weightsRaster.GetGeoTransform()
            geoBase = baseImg.GetGeoTransform()
            reproject = [1 for a, b in zip(geoSmooth, geoBase) if not a == b]
            if reproject:
                weightsRaster = rt.warpRaster(weightsRaster, baseImg,
                                            resampleOption='nearest',
                                            outputURI=None, outFormat='MEM')
        
        baseImg = None
        nodataWeights = weightsRaster.GetRasterBand(1).GetNoDataValue()
    
    # Prepare an empty dictionary to hold the results
    average = {}
    
    # Loop through the images to compute the average for each
    for img, date, mask in zip(images, datesImg, maskedQual):
        
        # Import the image
        baseImg = gdal.Open(img)
        
        if mask:
            maskImg = gdal.Open(mask)
        
        # Loop through the blocks to compute the average for each
        
        # Get the size of the rasters to identify the limits of the blocks to loop through
        band = baseImg.GetRasterBand(1)
        # Get the no data value
        nodataBase = band.GetNoDataValue()
        # Get the size of the raster
        xsize = band.XSize
        ysize = band.YSize
        # Get the number of blocks in x and y directions based on the block size
        xBlocks = int(round(xsize / blockXSize)) + 1
        yBlocks = int(round(ysize / blockYSize)) + 1
        band = None
        
        sumNdvi = []
        sumWeights = []
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                # Read the block from the image
                blockBase = rt.readRasterBlock(baseImg,
                                            xStep * blockXSize, yStep * blockYSize,
                                            blockXSize, blockYSize)
                
                # Recast the type to be sure
                blockBase = blockBase.astype(np.float32)
                
                # Replace the no data values by 0
                blockBase[np.logical_or(
                    np.logical_or(blockBase == nodataBase, np.isnan(blockBase)),
                    np.logical_or(blockBase > 1, blockBase < -1))] = 0.
                
                if weightsRaster:
                    # Read the block from the weights
                    blockWeight = rt.readRasterBlock(weightsRaster,
                                                  xStep * blockXSize, yStep * blockYSize,
                                                  blockXSize, blockYSize)
                    
                    # Recast the type to be sure
                    blockWeight = blockWeight.astype(np.float32)
                
                    # Replace the no data values by 0
                    blockWeight[np.logical_or(
                        np.logical_or(blockWeight == nodataWeights, np.isnan(blockWeight)),
                        np.logical_or(blockWeight > 1, blockWeight < 0))] = 0.
                    
                    # Set the no data values in the image as zero weight
                    blockWeight[blockBase == 0.] = 0.
                    
                else:
                    blockWeight = np.copy(blockBase)
                    
                    # Replace the non 0 data values by 1
                    blockWeight[blockWeight != 0.] = 1.
                
                # Remove the weights were the pixels are of low quality
                if mask:
                    maskBlock = rt.readRasterBlock(maskImg,
                                            xStep * blockXSize, yStep * blockYSize,
                                            blockXSize, blockYSize)
                    blockWeight[maskBlock == -3000] = 0.
                
                # Estimate the weighted sum for each pixel
                sumNdvi.append(np.sum(np.multiply(blockBase, blockWeight)))
                sumWeights.append(np.sum(blockWeight))
                
        # Combine for the entire image
        sumNdvi = np.sum(sumNdvi) / np.sum(sumWeights)
        
        # Add to the output dictionary along with the date
        average[date.strftime('%Y-%m-%d')] = sumNdvi
    
        # Close the raster
        baseImg = None
    
    return(average)


def maskQualityVI(ndviRaster, qualityRaster, outRaster=None, nodataOut=None):
    '''
    Function to mask a composite modis ndvi image with a quality layer.
    
    Each pixel of the quality layer is converted to 16 bit to screen
    quality problems and mask pixels with low quality.
    
    ndviRaster (str): Full address of the raster with the ndvi or other value
        of interest.
    qualityRaster (str): Full address of the raster with the quality 
        information as provided in MODIS data.
    outRaster (str): Full address for the output raster. If None, the input 
        raster is overwritten.
    nodataOut (num): value to use for data of low quality. Can be different 
        from the no data value of the raster. If None, the no data value of the
        raster is used.
    '''
    
    # Import the quality information image
    q = gdal.Open(qualityRaster)
    
    # Transform into array
    qArray = q.GetRasterBand(1).ReadAsArray()
    
    # Get the type
    qType = gdal.GetDataTypeName(q.GetRasterBand(1).DataType)
     
    # Import the modis image with the ndvi
    ndvi = gdal.Open(ndviRaster)
    # Transform into array
    ndviArray = ndvi.GetRasterBand(1).ReadAsArray()
    
    # Get the type
    ndviType = gdal.GetDataTypeName(ndvi.GetRasterBand(1).DataType)
    
    # Get the no data value
    nodataNdvi = ndvi.GetRasterBand(1).GetNoDataValue()
    nodataQual = q.GetRasterBand(1).GetNoDataValue()
    
    if not nodataOut:
        nodataOut = nodataNdvi
                    
    # Loop through the arrays to mask the ndvi
    
    for p, q in np.nditer([ndviArray, qArray],
                          flags=['reduce_ok'],
                          op_flags=[['readwrite'], ['readonly']],
                          op_dtypes=[ndviType, qType]):
        p[...] = screenQualityVI(p, q, nodataNdvi, nodataQual, nodataOut)
    
    if outRaster:
        try:
            os.remove(outRaster)
        except:
            pass
        
        # Create an empty raster copy
        out = rt.newRasterFromBase(ndvi, outRaster, 'GTiff', nodataNdvi,
                                   getattr(gdal,
                                           rt.getGDALTypeFromNumber(
                                               ndvi.GetRasterBand(1).DataType)))
        
        # Write to the new raster
        out.GetRasterBand(1).WriteArray(ndviArray)
        
        out.FlushCache()
        out = None
    else:
        ndvi.GetRasterBand(1).WriteArray(ndviArray)
        ndvi.FlushCache()
    
    # Close the rasters
    ndvi = None
    q = None
    
    
def screenQualityVI(pixel, index, nodataP, nodataI, nodataMask):
    '''
    Function to screen the quality of a pixel of ndvi in composite 
        modis product based on quality layer.
    
    pixel (num): value of the pixel of the raster to mask
    index (num): 16bit number of the quallity layer pixel
    nodataP (num): value of the no data for the raster to mask
    nodataI (num): value of the no data for the quality raster
    nodataMask (num): value to use for the pixel of low quality
    '''
    
    # If no quality information return the no data value of the input
    if index == nodataI or pixel == nodataP:
        return nodataP
    
    bit16 = np.binary_repr(index, width=16)
    
    if (bit16[2:5] in ['001', '010', '100', '101'] and bit16[14:16] in ['00', '01', '10'] and 
        bit16[9:13] in ['0000', '0001', '0010', '0011', '0100', '0101', '0111',
                        '1000']):
        p = pixel
    
    elif bit16[2:5] in ['000', '011', '110', '111']:
        p = nodataP
    else:
        p = nodataMask
    
    return(p)


def percMissingStack(images, outName, nodata=None):
    '''
    Takes a list of images/rasters of different dates and returns a 
    single raster wiwth for each pixel the % of missing
    
    A no data value should be specified if it is missing from the rasters
    The no data value should be the same for all the rasters
    
    images (list): List of full addressesof the images
    outName (str): Full address for the output raster (each pixel will hold 
        a % of non missing values
    nodata (num): no data value. Value of the missing data.
    '''
    
    # Import all the images
    toProcess = [gdal.Open(f) for f in images]
        
    # Get the no data value if any
    i = 0
    nodataImg = toProcess[0].GetRasterBand(1).GetNoDataValue()
    while not nodataImg and i < len(toProcess) - 1:
        i += 1
        nodataImg = toProcess[i].GetRasterBand(1).GetNoDataValue()
    
    if not nodataImg:
        nodataImg = nodata
    
    if not nodataImg:
        print('The nodata information is missing from all the images')
        return(False)
    
    # Create an empty raster to host the final percentages
    out = rt.newRasterFromBase(toProcess[0], outName, 'GTiff', np.nan, gdal.GDT_Float32)
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = toProcess[0].GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Set the block size
    blockXSize = 256 * 4
    blockYSize = 256 * 4
    # Get the number of blocks in x and y directions
    xBlocks = int(round(xsize / blockXSize)) + 1
    yBlocks = int(round(ysize / blockYSize)) + 1
    
    totBlocks = xBlocks * yBlocks
    progress = 1
        
    for xStep in range(xBlocks):
        for yStep in range(yBlocks):
            print('Processing block ' + str(progress) + ' of ' + str(totBlocks))
            progress += 1
            
            blocks = [rt.readRasterBlock(p, xStep * blockXSize, yStep * blockYSize, blockXSize, blockYSize) for p in toProcess]
            
            # Bring the blocks together into one single array
            blocks = np.dstack(blocks)
            
            # Do the smoothing for each pixel
            blocks = percMissing(blocks, nodataImg)
            
            # Change the values in the output raster
            out.GetRasterBand(1).WriteArray(blocks, xStep * blockXSize, yStep * blockYSize)
        
    # Close the rasters
    for p in toProcess:
        p = None
    out = None


def percMissing(block, nodata):
    '''
    Function to get the share of missing data.
    The function will loop through all the columns in the first two dimensions, 
    counting the missing in the third dimension.
    It returns an array of the share of missing.
    
    block (numpy array): Array to be processed
    nodata (int): value for the no data in the input array.
    '''
    
    extent = block.shape
    outBlock = np.ones(extent[0:2], dtype=float)
    tot = float(extent[2])
    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Check if entire pixel vector is nan, just return it
            if np.all(pixel == nodata) or np.all(np.isnan(pixel)):
                continue
            
            outBlock[X, Y] = (pixel == nodata).sum() / tot
    
    return outBlock

def plotModisLtavg(inDic, ltAvgStart, ltAvgEnd, dateStartChart, yearsPlot, outFolder):
    '''
    Function to compute the difference of the values in inDic to the 
    long term average (by date).
    
    inDic (dict): Dictionary of dictionaries of the data to plot. 
            The keys of the main dictionary are the dates and the values are 
            dictionaries. These secondary dictionaries have the variable to plot as 
            keys and the value for that date as values.
            The dates should be string, with the format 'yyyy-mm-dd'
            All dates should contain the same variables.
    ltAvgStart (num): start year for computing the long term average (included)
    ltAvgEnd (num): end year for hte long term average (included)
    dateStartChart (str): Date of the year at which the chart will start.
            Format 'mm-dd'
    yearsPlot (list): list of years to plot. Each year will be a line.
            If the start date for the chart is not the beginning of the year, 
            each year will include the beginning dates of the following year.
    outFolder (str): Full address of the output folder where to export the charts
    '''
    
    ##Take the input dictionary and extract the dates and the values
    
    #Get the dates in the input data
    inDates = inDic.keys()
    #Transform the dates into date format
    inDates = [datetime.strptime(d, '%Y-%m-%d').date() for d in inDates]
    #Sort the dates
    inDates.sort()
    
    #Create empty dictionary to host the lists of values
    #Export them in the same order as the dates
    values = {}
    for d in inDates:
        for k,v in inDic[d.strftime('%Y-%m-%d')].iteritems():
            #Add the new variable to the dictionary
            if not k in values:
                values[k] = []
            #Add the value
            values[k].append(v)
    
    #Check that all the dates contain the same variables
    lSeries = [len(v) for v in values.itervalues()]
    if not lSeries.count(lSeries[0]) == len(lSeries):
        print('Dates in input dictionary do not all contain the same variables')
        return False
    lSeries = None
    
    #Transform each date into a day of the year
    inDays = [int(d.strftime('%j')) for d in inDates]
    #Extract the year for each date
    inYears = [int(d.strftime('%Y')) for d in inDates]
    
    #Get the unique days
    uniqueDays = list(set(inDays))
    uniqueDays.sort()
    
    #Transform the start date for the charts into date format
    dateStartChart = datetime.strptime(dateStartChart, '%m-%d').date()
    #Transform into day of the year
    dateStartChart = int(dateStartChart.strftime('%j'))
    #Check if there are days in the dates that are less than the start date
    lessDays = any(d<dateStartChart for d in uniqueDays)
    
    #Import the color scheme
    if len(yearsPlot)<11:
        colormap = cm.tab10.colors
    elif len(yearsPlot)<21:
        colormap = cm.tab20.colors
    else:
        print('The code cannot map more than 20 variables in a chart')
        return False
    
    #Loop through each variable
    for V in values:
        #For each date compute the long term average
        ltavgs = {}
        for D in uniqueDays:
            #Get all the values for that day across all the years
            ltavg = [v for v,d,y in zip(values[V],inDays,inYears) 
                     if (d==D and y>=ltAvgStart and y<=ltAvgEnd)]
            #Compute the long term average for that day and add to list
            ltavgs[D] = sum(ltavg)/float(len(ltavg))
        
        #Compute the difference to the long term average for that value
        values[V] = [v-ltavgs[d] for v,d in zip(values[V],inDays)]
        
        #Create the lists for the plotting
        toPlot = {}
        #Save the labels and dictionary keys in order
        labelKeys = []
        #Loop through the years to plot
        for Y in yearsPlot:
            if lessDays:
                k = str(Y)[2:]+'/'+str(Y+1)[2:]
            else:
                k = str(Y)[2:]
            labelKeys.append(k)
        
            toPlot[k] = [
                [v for v,d,y in zip(values[V],inDays,inYears) 
                         if ((d>=dateStartChart and y==Y) or (d<dateStartChart and y==Y+1))],
                [d for d,y in zip(inDays,inYears) 
                         if ((d>=dateStartChart and y==Y) or (d<dateStartChart and y==Y+1))]
            ]
            
            #Transforms the days to strings for the plotting
            toPlot[k][1] = [datetime(2006+(d<dateStartChart), 1, 1) + timedelta(d - 1) for d in toPlot[k][1]]
            
        #Create the plot for that variable
        fig, ax = plt.subplots()
        fig.set_size_inches(10,7) #In inches... cm/2.54
        
        i = 0
        for p in labelKeys:
            ax.plot(toPlot[p][1], toPlot[p][0], color=colormap[i], label=p)
            i += 1
        #Add horizontal line at 0
        plt.axhline(color='k')
        #Place legend
        ax.legend(loc='center right', bbox_to_anchor=(1.09, 0.22), 
                  title='Data Year', facecolor='white', 
                  framealpha=1) #ncol=3
        #plt.legend(loc='upper left', title='Data Year')
        #Format the dates
        fig.autofmt_xdate()
        ax.xaxis.set_major_formatter(dates.DateFormatter("%b"))
        ax.xaxis.set_major_locator(dates.MonthLocator())
        #Add grid
        plt.grid()
        #Change color of background
        ax.set_facecolor('whitesmoke')
        #Add labels
        plt.xlabel('Date', weight='bold')
        plt.ylabel('Difference to Long Term Average CHI (>0 is better than average)', weight='bold')
        plt.title('Crop Health Index Annual Profiles (Diff. to Long Term Average): '+V, weight='bold')
        #Fit layout for smaller image
        plt.tight_layout()

        #Save to disk
        plt.savefig(outFolder+'/'+V+'_CHI_annual_profiles_diff_to_long_term.png', dpi=100)
        
def mapModisRanking(mapSize, mapTitle, mapFile,
                    boundaryFile, notePosition, 
                    noteSize, legendPosition, 
                    legendSizes, outName, 
                    outRes=200, backgroundLabel='',
                    citiesFile=None, citiesField=None,
                    citiesLabelSize=None, 
                    citiesMarkerSize=None, 
                    scaleLen=100, scaleSize=12, 
                    scalePosition=(0.6,0.01)):
    
    '''
    Maps the ranking raster and export to file
    
    mapSize (tuple of num): Size of the ap in inches
    mapTitle (str): Title of the map 
    mapFile (str): Full address of the raster with 
        the ranking information
    boundaryFile (str): Full address of the shapefile 
        with the boundaries information
    notePosition (tuple of num): Two values between 0 and 
        1 giving the position of the top left corner of the 
        note on the chart. 
        The value are relative to the chart, so (0,0) is the
        bottom left corner and (1,1) is the top right corner. 
    legendPosition (tuple of num): Two values giving the 
        position of the legend. Same as notePosition
    outName (str): Full address of the output name for the
        chart (.png)
    outRes (int): Dpi resolution for the output chart
    backgroundLabel (str): Label for the white background
        pixels
    citiesFile (str): Full address of the shapefile with the 
        information on the cities.
    citiesField (str): Name of the field in the shapefile with 
        name of the cities to use as label
    citiesLabelSize (int): Size of the labels
    '''

    #Create new figure window
    fig = plt.figure(figsize=mapSize)  # a new figure window
    ax = fig.add_subplot(1, 1, 1)  # specify (nrows, ncols, axnum)
    
    #Remove frame of subplot
    ax.axis('off')
    
    #Add title
    ax.set_title(mapTitle, fontsize=12, 
                 weight='bold', y=1.02)
    
    # Read the data and metadata
    datafile = gdal.Open(mapFile)
    bnd1 = datafile.GetRasterBand(1).ReadAsArray()
    
    #Get no data value
    nodata = datafile.GetRasterBand(1).GetNoDataValue()
    
    #Change data type and remove no data value from raster
    bnd1 = bnd1.astype(float)
    bnd1[bnd1==nodata] = np.nan
    
    #Get raster size, projection, and resolution
    nx = datafile.RasterXSize # Raster xsize
    ny = datafile.RasterYSize # Raster ysize
    gt = datafile.GetGeoTransform()
    xres = gt[1]
    yres = gt[5]
    
    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * nx) - xres * 0.5
    ymin = gt[3] + (yres * ny) + yres * 0.5
    ymax = gt[3] - yres * 0.5
    
    # create a grid of lat/lon coordinates in the original projection
    (lon_source,lat_source) = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
    
    #Create the basemap
    edge = ((xmax-xmin)*0.03,(ymax-ymin)*0.03) #Add an edge around the raster
    mapR = Basemap(projection='cyl',llcrnrlat=ymin-edge[1],urcrnrlat=ymax+edge[1],\
                llcrnrlon=xmin-edge[0],urcrnrlon=xmax+edge[0] , resolution='i', ax=ax)
    
    #Prepare the color map
    #Colors and values for scale
    colorsScale = [(255,254,141),(239,48,166),(197,38,182),(113,28,198),(19,0,236),(28,67,198),(5,251,255)]
    stopsScale = [-3,10,30,50,70,90,110]
    #Normalize the colors to 0-1
    norm = mpl.colors.Normalize(vmin=0.,vmax=255.)
    colorsScale = [tuple(norm(v) for v in T) for T in colorsScale]
    #Normalize the values for the scale
    norm = mpl.colors.Normalize(vmin=-3.,vmax=110.)
    stopsScale = [norm(v) for v in stopsScale]
    #Combine into color scale for mapping
    purpleBlue = zip(stopsScale,colorsScale)
    #Create the segmented color map
    purpleBlueLinear = col.LinearSegmentedColormap.from_list('purpleBlue',purpleBlue, N=256, gamma=1.0)
    #Create labels for colormap
    purpleBlueLabels = ['Below min', 
                        'Min of the 10 ref. years',
                        '',
                        'Median of the 10 ref. years',
                        '',
                        'Max of the 10 ref. Years',
                        'Above max']
    
    # project in the original Basemap and plot with pcolormesh
    mapR.pcolormesh(lon_source,lat_source,bnd1.T, cmap=purpleBlueLinear)
    
    #Add boundary
    
    aoi_info = mapR.readshapefile(re.sub('.shp$','',boundaryFile), 
                                  'aoi', color='black',
                                  linewidth=1.2)
    
    if citiesFile:
        #Add major cities
        cities_info = mapR.readshapefile(re.sub('.shp$','',citiesFile), 
                                         'cities')
        
        #Add the city names as labels
        cityFont = {'fontname':'Arial', 'size':str(citiesLabelSize), 
                    'color':'black', 'weight':'bold'}
        for info, city in zip(mapR.cities_info, mapR.cities):
            mapR.plot(city[0], city[1], marker='o', 
                      color='black', markersize=citiesMarkerSize, 
                      markeredgewidth=2) #'o' for circle, '.' for point
            plt.text(city[0]+0.015, city[1]+0.015, 
                     unicode(info[citiesField], 'utf-8'), 
                     path_effects=[pe.withStroke(linewidth=2, 
                                                 foreground="white")], 
                     **cityFont)
            #The path effect adds a buffer around the letters for legibility
    
    #Add scale bar
    scaleParam = {'startLon': xmin+(xmax-xmin)*scalePosition[0], 
                  'startLat': ymin+(ymax-ymin)*scalePosition[1],
                  'lengthKm': scaleLen, 
                  'yoffset': 0.02}
    #Get initial lon lat in map units
    lon1,lat1 = mapR(scaleParam['startLon'],scaleParam['startLat'],inverse=True)
    #Get final lon lat from distance
    gc = pyproj.Geod(a=mapR.rmajor,b=mapR.rminor)
    lon2, lat2, az = gc.fwd(lon1,lat1,90,scaleParam['lengthKm']*1000)
    #Get back the final lon lat in map units
    x2,y2 = mapR(lon2,lat2,inverse=False)
    #Plot the lines for the scale
    barHeight = abs(scaleParam['startLon']-x2)/100.
    mapR.plot([scaleParam['startLon'],x2],
              [scaleParam['startLat'],scaleParam['startLat']],color='k')
    mapR.plot([scaleParam['startLon'],scaleParam['startLon']],
              [scaleParam['startLat']-barHeight,
               scaleParam['startLat']+barHeight],color='k')
    mapR.plot([x2,x2],
              [scaleParam['startLat']-barHeight,
               scaleParam['startLat']+barHeight],color='k')
    scaleFont = {'fontname':'Arial', 'size':scaleSize, 'color':'black', 'weight':'normal', 
                 'horizontalalignment':'center'}
    plt.text(scaleParam['startLon'],scaleParam['startLat']+barHeight+0.01,
             '0', **scaleFont)
    plt.text(x2,scaleParam['startLat']+barHeight+0.01,
             '%s km' %(scaleParam['lengthKm']),    
             **scaleFont)
    
    
    #Add note
    commentFont = {'fontname':'Arial', 'size':noteSize, 'color':'black'}
    comment = 'Notes:'\
        '\n--Most pixels contain other land uses \nbesides coffee.'\
        '\n--The index shows health of vegetation in \neach pixel compared to reference years, \n'\
        'not coffee production directly. '\
        '\n--Vegetative health is affected mostly by \n'\
        'natural factors (rain, etc.) but can also be \naffected by human intervention (pruning, \netc.).'
    plt.text(notePosition[0], notePosition[1], comment, 
             horizontalalignment='left',
             verticalalignment='center',
             transform = fig.transFigure, #ax.transAxes for position relative to axes 
             **commentFont)
    
    #Prepare color legend
    if backgroundLabel:
        cmap = mpl.colors.ListedColormap([(1,1,1)]+colorsScale)
        #bounds = range(len(colorsScale)+1)
    else:
        cmap = mpl.colors.ListedColormap(colorsScale)
    bounds = range(cmap.N+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    #Add axis for the color legend
    legendFont = {'fontsize': legendSizes[1],
                  'fontweight': 'normal',
                  'verticalalignment': 'center'}
    axColors = fig.add_axes([legendPosition[0],legendPosition[1], 0.03, 0.18]) #left, bottom, width, height] in fractions of figure width and height
    cb = mpl.colorbar.ColorbarBase(axColors, cmap=cmap,
                                    norm=norm,
                                    boundaries=bounds,
                                    ticks=[y+0.5 for y in range(cmap.N+1)],
                                    spacing='uniform',
                                    orientation='vertical')
    cb.ax.set_title('Legend', fontsize=legendSizes[0], weight='bold', x=0.7, y=1.05)
    axColors.tick_params(axis=u'both', which=u'both',length=0)
    if backgroundLabel:
        legendLabels = [backgroundLabel]+purpleBlueLabels
    else:
        legendLabels = purpleBlueLabels
    axColors.set_yticklabels(legendLabels, fontdict=legendFont)
    
    #Fit layout for smaller image
    plt.tight_layout()
    
    #Save plot
    plt.savefig(outName,dpi=outRes)