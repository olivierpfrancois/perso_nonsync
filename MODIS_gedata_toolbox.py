# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Script to download MODIS images from the server, 
    and mosaic the images and mask them for each region
    to process.
 """

# Some commonly used imports

import pymodis as pm
import re, os, subprocess
from datetime import datetime, timedelta
from osgeo import gdal, gdalconst, ogr
import numpy as np
from csv import DictWriter

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# #FUNCTIONS


def downloadMODIS(dstFolder, pwd, user, tiles, product, startDownload=None, endDownload=None):
    '''
    Function to download modis images from server. Returns a list of the names of the newly downloaded files.
    
    dstFolder (str): Address of the destination folder where to save the downloaded images
    pwd (str): Password for the e4ftl01.cr.usgs.gov website
    user (str): Username for the e4ftl01.cr.usgs.gov website
    tiles (list of str): List of  tiles to download (e.g. 'h13v10')
    product (str): Product to download from server (e.g. 'MOD13Q1.006')
    startDownload (str): Start date for the product download (format YYYY-MM-DD). 
                        If None, will default to date of most recent MODIS file on disk if any, or stop the process
    endDownload (str): End date for the product download (format YYYY-MM-DD). If None, defaults to today
    '''
    
    if not endDownload:
        # Defaults to today's date
        endDownload = datetime.now()
        endDownload = endDownload.strftime('%Y-%m-%d')
        
    # Get the names of the hdf files already downloaded and processed
    # thereTif = [f for f in os.listdir(dstFolder) if f.endswith('.tif')]
    thereHdf = [f for f in os.listdir(dstFolder) if f.endswith('.hdf')]
    # Dates of these files
    if thereHdf:
        datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in thereHdf]
        datesHdfD = [datetime.strptime('0101' + d[0:4], "%d%m%Y").date() + timedelta(days=int(d[4:])) for d in datesHdf]
        # datesTif = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})_', f).group(1) for f in thereTif]
        # datesTif = [datetime.strptime(d, "%Y-%m-%d").date() for d in datesTif]
    else:
        datesHdf = []
        # datesTif = []
    
    # Get latest date of the file on disk
    if datesHdf and not startDownload:
        startDownload = max(datesHdfD).strftime('%Y-%m-%d')
    
    if not startDownload:
        print('There are no files on disk. Choose a start date for the download.')
        return
    
    # Download 
    down = pm.downmodis.downModis(destinationFolder=dstFolder, password=pwd, user=user,
                url="https://e4ftl01.cr.usgs.gov", tiles=tiles, path='MOLT', product=product,
                today=endDownload, enddate=startDownload, jpg=False, debug=False, timeout=30, checkgdal=True)
    down._connectHTTP()
    down.downloadsAllDay()
    
    # Get the list of new hdf files on disk
    newHdf = [f for f in os.listdir(dstFolder) if f.endswith('.hdf')]
    newHdf = [f for f in newHdf if f not in thereHdf]
    
    if newHdf:
        print(str(len(newHdf)) + ' images downloaded')
        
        # Remove the .txt files
        xml = [f for f in dstFolder if f.endswith('.txt') or f.endswith('.log')]
        for f in xml:
            os.remove(os.path.join(dstFolder, f))
    else:
        print('No new images downloaded')
    
    return newHdf

   
def mosaicMODISWrapper(root, srcFolder, tmpFolder, regions, regionsOut, regionsBoundaries, tiles, subset, suffix, startMosaic=None, endMosaic=None):
    '''
    Function to mosaic the tiles of modis images and clip them to a series of regions.
    Does not return anything.
    
    root (str): Address of root folder where  regions folders are located
    srcFolder (str): Full address where the input downloaded tiles are located
    tmpFolder (str): Full address of folder where the temporary mosaic before cutting to regions will be stored. 
            !!!!! THIS FOLDER SHOULD BE EMPTY TO START, IT WILL BE EMPTIED AT THE END
    regions (list of str): Names of the regions to process. Each region should have a folder inside root with the same name.  
    regionsOut (str): Name of the folder inside the regions folders where the mosaiced and clipped images should be stored. 
            It should be the same for all the regions
    regionsBoundaries (list of str): List the full address of the shapefiles with the regions boundaries
    tiles (list of str): List of  tiles to mosaic (e.g. 'h13v10')
    subset (str): string with the layers to extract from the hdf images (e.g. '1 0 0 0 0 0 0 0 0 0 0 0')
    suffix (list): list of string, suffixes to assign for naming each of the subsets
    startMosaic (str): Starting date for the files to mosaic. If None, will process all the files found on the disk.
    endMosaic (str): Ending date for the files to mosaic. If None, defaults to today
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
    datesHdfD = [datetime.strptime('0101' + d[0:4], "%d%m%Y").date() + timedelta(days=int(d[4:]) - 1) for d in datesHdf]
    
    # Get the dates to process
    if startMosaic:
        dates = {d for d, D in zip(datesHdf, datesHdfD) if D >= startMosaic and D <= endMosaic}
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
        # dlong = datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])-1)
        
        # Progress
        print('Processing date ' + dlong.strftime('%Y-%m-%d'))
        
        # Get the file names for that date
        files = [f for f in thereHdf if d in f]
        
        # Check if all the tiles were completed
        complete = True
        for t in tiles:
            if not complete:
                break
            
            complete = any(t in f for f in files)
        
        if not complete:
            print('Not all tiles were downloaded for date ' + dlong.strftime('%Y-%m-%d') + ". Cannot process.")
            continue
        
        # Create the output name for the mosaic file
        outName = re.sub('\..+', '', files[0])  # Get the product name
        outName = outName + '_' + dlong.strftime('%Y-%m-%d') + '_h' + '-'.join(map(str, tilesH)) + 'v' + '-'.join(map(str, tilesV)) + '_250m_16_days'
        
        mosaicMODIS(images=[os.path.join(srcFolder, f) for f in files],
                    subset=subset,
                    suffixes=suffix,
                    tempFolder=tmpFolder,
                    outFile='/'.join([tmpFolder, outName + '.tif']))
        
        # Loop through the regions to clip and mask the mosaics
        for s, shp in zip(regions, regionsBoundaries):
            
            # Loop through the mosaics
            for l in range(len(suffix)):
                
                # Clip the mosaic by the extent of the the shapefile
                clipMaskRasterByShp(shp=shp,
                                    raster='/'.join([tmpFolder, outName + '_' + suffix[l] + '.tif']),
                                    outRaster=root + '/' + s + '/' + regionsOut + '/' + outName + '_' + suffix[l] + '.tif',
                                    nodata=-3000)
                
                # Remove intermediary file
                os.remove('/'.join([tmpFolder, outName + '_' + suffix[l] + '.tif']))
        
        # Remove the .xml, .hdf and .vrt intermediary files
        xml = [f for f in os.listdir(os.path.join(tmpFolder)) if f.endswith('.xml') or 
                    f.endswith('.hdf') or f.endswith('.vrt')]
        for f in xml:
            os.remove(os.path.join(tmpFolder, f))


def mosaicMODIS(images, subset, suffixes, tempFolder, outFile):
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


def clipMaskRasterByShp(shp, raster, outRaster, nodata=None):
    
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
    ulX, ulY = world2Pixel(geoTrans, minX, maxY)
    lrX, lrY = world2Pixel(geoTrans, maxX, minY)
    
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
    
    nodataImg = band.GetNoDataValue()
    if not nodataImg:
        nodataImg = nodata
    
    if not nodataImg:
        print('No NODATA information in raster and none provided')
        return False
    
    # Extract the image information in that extent
    # using readAsArray xoffset, yoffset, xextent, yextent
    clip = band.ReadAsArray(ulX, ulY, int(lrX - ulX), int(lrY - ulY)).astype(np.float)

    # Remove rasters to save
    if os.path.isfile(outRaster):
        os.remove(outRaster)
    
    # Create an empty raster with the right size
    rows, cols = clip.shape
    driver = gdal.GetDriverByName("GTiff")
    new_raster = driver.Create(outRaster, cols, rows, 1, getattr(gdal, getGDALTypeFromNumber(dType)),
                               options=['COMPRESS=LZW', 'BIGTIFF=IF_NEEDED'])
    new_raster.SetProjection(projection)
    new_raster.SetGeoTransform(geoRegion)

    new_raster.GetRasterBand(1).SetNoDataValue(nodataImg)
    new_raster.GetRasterBand(1).Fill(nodataImg)

    # Rasterize the region into that raster
    gdal.RasterizeLayer(new_raster, [1], maskLayer, burn_values=[1.], options=['ALL_TOUCHED=FALSE'])
    
    # Get the resulting raster band as an array
    new = new_raster.GetRasterBand(1).ReadAsArray().astype(np.float)
    new[new == 1.] = clip[new == 1.]
    
    # Export the values
    new_raster.GetRasterBand(1).WriteArray(new, 0, 0)
    
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
    
    image = None


def smoothMODIS(root, regions, regionsIn, regionsOut, startSmooth, endSmooth, regWindow, avgWindow,
                startSaveSmooth=None, endSaveSmooth=None, algo='Swets'):
    '''
    Function to do a temporal smoothing of a time series of identical images. Does not return anything, saves smoothed images on disk.
    The input files should contain the date in their name in the format '_([0-9]{4}-[0-9]{2}-[0-9]{2})'
    The processed files will have the same name plus a 'smooth_' prefix
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region should have a folder inside root with the same name.
    regionsIn (str): Name of the folder inside the regions folders where to find the input files to be smoothed 
            It should be the same for all the regions 
    regionsOut (str): Name of the folder inside the regions folders where to save the smoothed images 
            It should be the same for all the regions
    startSmooth (str): Starting date for the files to mosaic. If None, will process all the files found on the disk.
    endSmooth (str): Ending date for the files to mosaic. If None, defaults to today
    regWindow (int): size of the regression window (see Swets et al. for details)
    avgWindow (int): sie of the averaging window (see Swets et al. for details)
    startSaveSmooth (str): Starting date for the files to save. If None, will save all the processed files
    endSaveSmooth (str): Ending date for the files to save. If None, will save all the processed files
    algo(str): Name of algorithm to use for the smoothing. Can be Swets or Savitzky
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
    
    # Loop through the regions to do the smoothing for each
    for r in regions:
        print('Processing region ' + str(r) + '...')
        
        # Import all the raw modis images on disk
        onDisk = [f for f in os.listdir(os.path.join(root, r, regionsIn)) if f.endswith('.tif')]
        
        # Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        # Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        
        # Keep only the files and dates within the dates to process
        onDisk = [f for f, d in zip(onDisk, datesAll) if d >= startSmooth and d <= endSmooth]
        datesAll = [d for d in datesAll if d >= startSmooth and d <= endSmooth]
        
        # Sort the two list by date
        datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
        # Create a mask to know which of the dates to save
        toSave = [d >= startSaveSmooth and d <= endSaveSmooth for d in datesAll]
        
        # Smooth the images
        smoothSeries(inRasters=onDisk, toSave=toSave,
                     outFolder=root + '/' + r + '/' + regionsOut,
                     algo='Swets', blockXSize=256, blockYSize=256)


def smoothSeries(inRasters, toSave, outFolder, algo='Swets', blockXSize=256, blockYSize=256):
    if not len(inRasters) == len(toSave):
        print('toSave should be a list of bollean with the same length as the inputs')
        return False
    
    # Import all the images to process in the right order
    toProcess = [gdal.Open(f) for f in inRasters]
    
    # Get the no data value if any
    nodata = toProcess[0].GetRasterBand(1).GetNoDataValue()
    
    # Remove rasters to save if any
    for s, f in zip(toSave, onDisk):
        if s and os.path.isfile(outFolder + '/smooth_' + os.path.basename(f)):
            os.remove(outFolder + '/smooth_' + os.path.basename(f))
    
    # Create empty copies of these files to use for the smoothed data only for the files to save
    processed = [new_raster_from_base(p, outFolder + '/smooth_' + os.path.basename(f),
                                      'GTiff', np.nan, gdal.GDT_Float32) 
                    if s else False for s, p, f in zip(toSave, toProcess, onDisk)]
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = toProcess[1].GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Get the number of blocks in x and y directions based on block size
    xBlocks = int(round(xsize / BlockXSize)) + 1
    yBlocks = int(round(ysize / BlockYSize)) + 1
    
    totBlocks = xBlocks * yBlocks
    progress = 1
    
    for xStep in range(xBlocks):
        for yStep in range(yBlocks):
            print('   ' + str(r) + ': Processing block ' + str(progress) + ' of ' + str(totBlocks))
            progress += 1
            
            blocks = [readRasterBlock(p, xStep * BlockXSize, yStep * BlockYSize, BlockXSize, BlockYSize) for p in toProcess]
            
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
            blocks = np.divide(blocks, 10000.)
            
            # Change the values in the output raster
            blocks = np.dsplit(blocks, len(toProcess))
            for s, p, b in zip(toSave, processed, blocks):
                if s:
                    p.GetRasterBand(1).WriteArray(b[:, :, 0], xStep * BlockXSize, yStep * BlockYSize)
    
    # Close the rasters
    for s, p in zip(toSave, processed):
        if s:
            p.FlushCache()
            p = None
    
    for p in zip(toSave, toProcess):            
        p = None


def smoothingSwets(block, regWindow, avgWindow, nodata=None):
    '''
    Function to do the temporal smoothing of a 3D numpy array.
    The function will loop through all the columns in the first two dimensions, smoothing the data in the 3rd dimension.
    The method is that of Swets et al. (1999, A weighted least-squares approach to temporal NDVI smoothing. ASPRS Annual Conference).
    It returns an array of the same shape with smoothed values, where the no data value is np.nan.
    
    The Swets et al. method runs a moving (weighted) regression window along the values in the input pixels. It then
    averages the predicted values using a second moving window the get the final smoothed values.
    
    block (numpy array): Array to be smoothed
    regWindow (int): size of the regression window (see Swets et al. for details)
    avgWindow (int): sie of the averaging window (see Swets et al. for details)
    nodata (int): value for the no data in the input array if any. If None, assumed to be np.nan.
            The output array will be returned with np.nan as the no data value.
    '''
    extent = block.shape
    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Check if entire pixel vector is nan, just return it
            if nodata and np.all(pixel == nodata):
                pixel.fill(np.nan)
                block[X, Y, :] = pixel
                continue
            elif np.all(np.isnan(pixel)):
                continue
            
            # Get the shape of the original data
            originShape = pixel.shape
            
            # Reshape the data
            pixel = np.reshape(pixel, (len(pixel), 1))
            pixel = pixel.astype(np.float32)
            
            # Replace no data values by nan if needed
            if nodata:
                pixel[pixel == nodata] = np.nan
            
            # Interpolate missing values if any
            # THIS SHOULD BE DONE BEFOREHAND FOR THE ENTIRE RASTER USING TEMPORAL AND SPATIAL INFERENCE
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
                # Local high values get a weight of 1.5, local middle values get a weight of 0.5 and local low values only 0.005
                if pixel[i - 1, 0] < pixel[i, 0] and pixel[i, 0] > pixel[i + 1, 0]:
                    weights[i, 0] = 1.5
                elif ((pixel[i - 1, 0] <= pixel[i, 0] and pixel[i, 0] <= pixel[i + 1, 0]) or (pixel[i - 1, 0] >= pixel[i, 0] and pixel[i, 0] >= pixel[i + 1, 0])):
                    weights[i, 0] = 0.5
                elif pixel[i - 1, 0] > pixel[i, 0] and pixel[i, 0] < pixel[i + 1, 0]:
                    weights[i, 0] = 0.005
            # For the last point
            if pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 2, 0] and pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 3, 0]:
                # If the last data point is greater than the previous 2, then assume it is a high point
                weights[len(pixel) - 1, 0] = 1.5
            elif pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 2, 0] and pixel[len(pixel) - 1, 0] < pixel[len(pixel) - 3, 0]:
                # If the last data point is greater than the previous only, then assume it is a middle point
                weights[len(pixel) - 1, 0] = 0.5
            elif (pixel[len(pixel) - 1, 0] < pixel[len(pixel) - 2, 0] and pixel[len(pixel) - 1, 0] >= pixel[len(pixel) - 3, 0] and 
                  weights[len(pixel) - 3, 0] != 0.005):
                # If less than the previous but more than the one before that is not a low point
                weights[len(pixel) - 1, 0] = 0.5
            else:
                # If the last point is less than the last 2, or greater than only one of the two, assume it is a low point
                weights[len(pixel) - 1, 0] = 0.002
            
            # Create a matrix with the data for this pixel and for the weights
            # For the data:
            # Each column will be the same
            dataRaw = np.repeat(pixel, len(pixel) + regWindow, axis=1)
            # Set to nan all data for each column besides the data used for the regressions
            ltri = np.tril_indices(n=len(pixel), k=-1, m=len(pixel) + regWindow)  # Lower triangle indices below the diagonal
            dataRaw[ltri] = np.nan
            utri = np.triu_indices(n=len(pixel), k=0, m=len(pixel))
            dataRaw[:, (regWindow):][utri] = np.nan
            # Remove the first two and last 3 columns, since they don't have enough points for a regression
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
                out = np.divide(np.absolute(y - np.mean(y)), np.std(y))  # Will be an outlier if greater than 3 (only extreme outliers should be picked up)
                out = np.reshape(out, y.shape)
                # Remove outliers before regression
                yout = y[out < 3]
                wout = w[out < 3]
                xout = x[out < 3]
                # Compute parameters of regression
                numerator = np.sum(wout) * np.sum(np.multiply(np.multiply(wout, xout), yout)) - np.sum(np.multiply(wout, yout)) * np.sum(np.multiply(wout, xout))
                denominator = np.sum(wout) * np.sum(np.multiply(wout, np.square(xout))) - np.square(np.sum(np.multiply(wout, xout)))
                b = np.divide(numerator, denominator)
                numerator = (np.sum(np.multiply(wout, yout)) - b * np.sum(np.multiply(wout, xout)))
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
                    res = np.mean(x[(int(np.floor(regWindow / 2.)) - t - 1):(int(np.floor(regWindow / 2.)) + t + 1)])
                else:
                    res = np.mean(x[(int(np.floor(regWindow / 2.)) - t):(int(np.floor(regWindow / 2.)) + t + 1)])    
                smoothed.append(res)
            
            smoothed = np.asarray(smoothed)
            smoothed = np.reshape(smoothed, originShape)
            
            block[X, Y, :] = smoothed
    
    return block


def smoothingSavitzky(block, nodata=None):
    '''
    Function to do the temporal smoothing of a 3D numpy array.
    The function will loop through all the columns in the first two dimensions, smoothing the data in the 3rd dimension.
    The method is that of Savitzky-Golay, using a window of 7 and a 2nd order polynomial.
    It returns an array of the same shape with smoothed values, where the no data value is np.nan.
    
    The Savitzky-Golay method runs a moving (weighted) window along the values in the input pixels. It simply 
    averages the values in the window to predict the output value.
    A window of 7 is still used for edages, but different weights are applied.
    
    block (numpy array): Array to be smoothed
    nodata (int): value for the no data in the input array if any. If None, assumed to be np.nan.
            The output array will be returned with np.nan as the no data value.
    '''
    extent = block.shape
    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Check if entire pixel vector is nan, just return it
            if nodata and np.all(pixel == nodata):
                pixel.fill(np.nan)
                block[X, Y, :] = pixel
                continue
            elif np.all(np.isnan(pixel)):
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
                # Replace the values that are lower in the original pixel by the values in the original smoothed pixel
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
            
            # Replace the smoothed pixel in the block
            block[X, Y, :] = pixel
    
    return block


def savitzkyAvg(pixel):
    # This function averages a pixel using a window 7 for the Savitzky Golay and returns a pixel of the same shape
    
    out = np.copy(pixel)
    
    for i in range(len(pixel)):
        if i < 3:
            # Take the moving window
            avg = pixel[0:7]
            # Do the weighted sum
            if i == 0:
                out[i] = np.sum(np.multiply(avg, np.array([0.1190, -0.0714, -0.1429, -0.0952, 0.0714, 0.3571, 0.7619])))
            elif i == 1:
                out[i] = np.sum(np.multiply(avg, np.array([-0.0714, -0.0000, 0.0714, 0.1429, 0.2143, 0.2857, 0.3571])))
            elif i == 2:
                out[i] = np.sum(np.multiply(avg, np.array([-0.1429, 0.0714, 0.2143, 0.2857, 0.2857, 0.2143, 0.0714])))
        elif i > (len(pixel) - 4):
            # Take the moving window
            avg = pixel[(len(pixel) - 7):]
            # Do the weighted sum
            if i == len(pixel) - 3:
                out[i] = np.sum(np.multiply(avg, np.array([0.0714, 0.2143, 0.2857, 0.2857, 0.2143, 0.0714, -0.1429])))
            elif i == len(pixel) - 2:
                out[i] = np.sum(np.multiply(avg, np.array([0.3571, 0.2857, 0.2143, 0.1429, 0.0714, -0.0000, -0.0714])))
            elif i == len(pixel) - 1:
                out[i] = np.sum(np.multiply(avg, np.array([0.7619, 0.3571, 0.0714, -0.0952, -0.1429, -0.0714, 0.1190])))
        else:
            # Take the moving window
            avg = pixel[i - 3:i + 4]
            # Do the weighted sum
            out[i] = np.sum(np.multiply(avg, np.array([-0.0952, 0.1429, 0.2857, 0.3333, 0.2857, 0.1429, -0.0952])))                
        
    return(out)
            

def createBaseline(root, regions, varieties, regionsIn, regionsOut, startRef, endRef, mask=None, outModelRaster=None):
    '''
    Function to create baseline images with the decile values over the period of interest. 
    Does not return anything, saves baselines images on disk (10 layer rasters).
    The input files should come from the smoothing process of modis data.
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region should have a folder inside root with the same name.
    varieties (list of lists): varieties to consider for each region. The function will create a series of baselines for 
            each regions and each varieties for that region. The name of the variety will be used to name the baselines.
            Elements in the list of lists are string, for example [['arabica','robusta'],['arabica'],...]
    regionsIn (str): Name of the folder inside the regions folders where to find the input modis rasters to 
            be used for the baselines. These files are outputs of the smoothing. 
            It should be the same for all the regions 
    regionsOut (str): Name of the folder inside the regions folders where to save the baselines images 
            It should be the same for all the regions
    startRef (str): Starting year for the files to consider for the baseline (included).
    endRef (str): Ending year for the files to consider for the baseline (included).
    mask (list of lists): None or list of lists of full addresses of the masks (.tif) to use for each of the regions and varieties.
            Each list element of the overall list should have the same length as the corresponding list for varieties or be None/empty.
            The masks should be at the same resolution and stack perfectly with the modis images.
            If None, no masking will be done.
    outModelRaster (list of lists): None or list of lists of full addresses of the rasters (.tif) to use as model for the output.
            or each of the regions and varieties.
        
    '''
    
    # Loop through the regions
    for r in range(len(regions)):
        # Get the names of all the smoothed rasters on file
        onDisk = [f for f in os.listdir(os.path.join(root, regions[r], regionsIn)) if f.endswith('.tif')]
        
        # Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        # Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        
        # Transform into days from start of the year and keep only the unique values
        days = {d.strftime('%j') for d in datesAll}
        # Transform back into a list and order by date
        days = [int(d) for d in days]
        days.sort()
        
        for v in range(len(varieties[r])):
            # Loop through the dates to create a baseline raster for each
            for d in days:
                # Get the names of all the rasters for this date
                dates = [datetime(y, 1, 1).date() + timedelta(days=d - 1) for y in range(startRef, endRef + 1, 1)]
                # dates = [datetime.strptime('0101'+str(y), '%d%m%Y').date()+timedelta(days=d-1) for y in range(startRef,endRef+1,1)]
                files = [f for f, date in zip(onDisk, datesAll) if date in dates]
                
                # Check if there is a mask and a model for the output raster
                if mask and mask[r] and mask[r][v]:
                    maskD = mask[r][v]
                else:
                    maskD = None
                if outModelRaster and outModelRaster[r] and outModelRaster[r][v]:
                    outModelRasterD = outModelRaster[r][v]
                else:
                    outModelRasterD = None
                
                # Prepare outName
                outName = (root + '/' + regions[r] + '/' + 
                           regionsOut + '/' + 'ndvi_deciles_0to100pct_ref_period-' + 
                           str(startRef) + '-' + str(endRef) + 
                           '_day' + str(d) + '_date-' + 
                           dates[0].strftime('%b-%d') + '_' + 
                           varieties[r][v] + '.tif')
                
                # Create decile raster
                createDecileRaster(images=files,
                              outFile=outName,
                              mask=maskD,
                              outModelRaster=outModelRasterD,
                              blockXSize=256, blockYSize=256)


def createDecileRaster(images, outFile, mask=None, outModelRaster=None, blockXSize=256, blockYSize=256):
    
    # Import all the images to use for estimating the deciles
    toProcess = [gdal.Open(os.path.join(root, regions[r], regionsIn, f)) for f in images]
    
    # Mask if there is a mask
    if mask:
        # I need to mask the rasters first to only have the coffee pixels when I change the resolution
        # Import the mask raster as an array
        p = gdal.Open(mask)
        nanMask = p.GetRasterBand(1).ReadAsArray()
        # Transform all the non zero values to nan
        nanMask = np.logical_or(nanMask == 0, np.isnan(nanMask))
        
        # Mask the rasters
        for p in toProcess:
            pBand = p.GetRasterBand(1).ReadAsArray()
            pBand[nanMask] = np.nan
            toProcess[-1].GetRasterBand(1).WriteArray(pBand)
        
        p = None
        pBand = None
        
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
        toProcess = [warp_raster(p, outModel, resampleOption='average', outputURI=None, outFormat='MEM') 
                     for p in toProcess]
    
    # Get the no data value if any
    nodata = toProcess[0].GetRasterBand(1).GetNoDataValue()
    
    # Remove existing raster if any
    if os.path.isfile(outFile):
        os.remove(outFile)
    
    # Create an empty copy with 10 layers to use for storing the deciles
    processed = new_raster_from_base(toProcess[0], outFile, 'GTiff',
                                     np.nan, gdal.GDT_Float32, bands=10)
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = toProcess[1].GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Get the number of blocks in x and y directions based on block size
    xBlocks = int(round(xsize / BlockXSize)) + 1
    yBlocks = int(round(ysize / BlockYSize)) + 1
    
    for xStep in range(xBlocks):
        for yStep in range(yBlocks):
        
            block = [readRasterBlock(p, xStep * BlockXSize, yStep * BlockYSize, BlockXSize, BlockYSize) for p in toProcess]
            
            # Bring the blocks together into one single array
            block = np.dstack(block)
            # Recast the type
            block = block.astype(np.float32)
            
            # Estimate the deciles for each pixel
            deciles = estimateDeciles(block, nodata)
            
            # Change the values in the output raster
            deciles = np.dsplit(deciles, 10)
            for i in range(10):
                processed.GetRasterBand(i + 1).WriteArray(deciles[i][:, :, 0], xStep * BlockXSize, yStep * BlockYSize)
    
    # Close the rasters
    for p in toProcess:
        p.FlushCache()
        p = None
    processed.FlushCache()
    processed = None


def estimateDeciles(block, nodata):
    extent = block.shape
    deciles = np.empty((extent[0], extent[1], 10))

    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X, Y, :])
            
            # Replace the no data value by nan if nodata was provided
            if nodata:
                pixel[pixel == nodata] = np.nan
            
            # Get the number of nan values in the pixel
            nbNodata = np.sum(np.isnan(pixel))
            
            # Return nan if there are not more than 6 valid values
            if len(pixel) - nbNodata < 6:
                deciles[X, Y, :].fill(np.nan)
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
    Function to rank each pixel of an image based on the baseline for the same day of year.
    Does not return anything, saves ranked image on disk.
    The input files should come from the smoothing and baseline functions.
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region should have a folder inside root with the same name.
    varieties (list of lists): varieties to consider for each region. The function will create a series of baselines for 
            each regions and each varieties for that region. The name of the variety will be used to name the baselines.
            Elements in the list of lists are string, for example [['arabica','robusta'],['arabica'],...]
    regionsIn (str): Name of the folder inside the regions folders where to find the input modis rasters to 
            be used for the ranking. These files are outputs of the smoothing. 
            It should be the same for all the regions 
    refDecilesIn (str): Name of the folder inside the regions folders where to find the input baseline rasters to 
            be used for the ranking. These files are outputs of the baseline creation. 
            It should be the same for all the regions 
    regionsOut (str): Name of the folder inside the regions folders where to save the baselines images 
            It should be the same for all the regions
    startRank (str): Starting date for the files to rank (included).
    endRank (str): Ending date for the files to rank (included).
    mask (list of lists): None or list of lists of full paths of the masks (.tif) to use for each of the regions and varieties.
            Each list element of the overall list should have the same length as the corresponding list for varieties or be None/empty.
            The masks should be at the same resolution and stack perfectly with the modis images.
            The mask should contain for each pixel either nodata or the density of the crop masked in the pixel (0 < d < 1).
            If None, no masking will be done.
    minDensity (decimal): None or list of values between 0 and 1. Minimum density of crop of interest to filter out pixels in the ranking. 
            All pixels with less than these densities will be masked, with one output per density.
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
        datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
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
                baseFile = [f for f in baseFiles if 'day' + str(day) + '_' and varieties[r][v] in f]
                baseFile = os.path.join(root, regions[r], refDecilesIn, baseFile[0])
                
                # Create the output name
                outName = os.path.isfile(root + '/' + regions[r] + '/' + 
                                         'ndvi_' + date.strftime('%Y-%m-%d') + 
                                         '_CompareToDecile_0BelowMin_110AboveMax_' + 
                                         varieties[r][v] + '.tif')
                
                # Rank the pixels of the image
                estimateRankRaster(image=os.path.join(root, regions[r], regionsIn, fileS),
                                   deciles=baseFile,
                                   densityMask=maskD,
                                   outFile=outName,
                                   minDensity=minDensity,
                                   BlockXSize=256, BlockYSize=256)


def estimateRankRaster(image, deciles, densityMask, outFile,
                       minDensity=None, BlockXSize=256, BlockYSize=256):
    
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
        smoothImg = warp_raster(smoothImg, baseImg, resampleOption='average', outputURI=None, outFormat='MEM')
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = baseImg.GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Get the number of blocks in x and y directions based on the block size
    xBlocks = int(round(xsize / BlockXSize)) + 1
    yBlocks = int(round(ysize / BlockYSize)) + 1
    band = None
    
    if not minDensity:
        minDensity = [0]
    
    for m in minDensity:
        if m >= 1 or m < 0:
            continue
        
        # Prepare output name
        outFile.replace('.tif', '_maskedbelow' + str(int(m * 100)) + '%.tif')
        
        # Remove existing raster if any
        if os.path.isfile(outFile):
            os.remove(outFile)
        # Create an empty copy for storing the deciles comparisons
        processed = new_raster_from_base(baseImg, outFile, 'GTiff',
                                         - 32768, gdal.GDT_Int16, bands=1)
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                # Read the block from the images
                blockSmooth = readRasterBlock(smoothImg, xStep * BlockXSize,
                                              yStep * BlockYSize,
                                              BlockXSize, BlockYSize)
                
                blockBase = [readRasterBlock(baseImg, xStep * BlockXSize,
                                             yStep * BlockYSize,
                                             BlockXSize, BlockYSize, band=b + 1) 
                             for b in range(baseImg.RasterCount)]
                
                # Read the block from the mask 
                blockMask = readRasterBlock(nanMask, xStep * BlockXSize,
                                            yStep * BlockYSize,
                                            BlockXSize, BlockYSize)
                
                if m > 0:
                    blockMask = blockMask < m
                
                # Bring the blocks together into one single array
                blockBase = np.dstack(blockBase)
                
                # Apply the mask
                blockSmooth[blockMask] = np.nan
                blockBase[blockMask] = np.nan
                
                # Recast the type to be sure
                blockSmooth = blockSmooth.astype(np.float32)
                blockBase = blockBase.astype(np.float32)
            
                # Estimate the placement for each pixel
                ranks = estimateRank(block=blockSmooth, ref=blockBase, nodata=nodata)
                
                # Change the values in the output raster
                processed.GetRasterBand(1).WriteArray(ranks, xStep * BlockXSize, yStep * BlockYSize)
    
    # Close the rasters
    smoothImg.FlushCache()
    smoothImg = None
    baseImg.FlushCache()
    baseImg = None
    processed.FlushCache()
    processed = None


def estimateRank(block, ref, nodata):
    extent = block.shape
    ranking = np.empty((extent[0], extent[1]), dtype=int)
    
    for X in range(extent[0]):
        for Y in range(extent[1]):
            testPixel = np.copy(block[X, Y])
            refPixel = np.copy(ref[X, Y, :])
            refPixel = np.reshape(refPixel, (len(refPixel), 1))
            
            # Replace the no data value by nan if nodata was provided
            if nodata:
                refPixel[np.logical_or(refPixel == nodata, np.isinf(refPixel))] = np.nan
                if testPixel == nodata or np.isinf(testPixel):
                    testPixel = np.nan
            else:
                refPixel[np.isinf(refPixel)] = np.nan
                if np.isinf(testPixel):
                    testPixel = np.nan
            
            # return nan if the test pixel is nodata
            if np.isnan(testPixel) or np.all(np.isnan(refPixel)):
                ranking[X, Y] = -32768
                continue
            
            # Get the position of the pixel in the reference deciles
            refPixel = np.reshape(refPixel, len(refPixel))
            
            rank = np.searchsorted(refPixel, testPixel) * 10
            
            # If the value is above the deciles, returns length of the vector
            if len(refPixel) == 10 and rank == 100:
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


def computeAvgNdvi(root, regions, varieties, regionsIn, avgWeights, weightField=None, startAvg=None, endAvg=None, alltouch=False):
    '''
    Function to compute the average ndvi value per region of interest.
    Does not return anything, saves averages in a txt file on disk.
    The input files should come from the smoothing function.
    
    root (str): Address of root folder where regions folders are located
    regions (list of str): Names of the regions to process. Each region should have a folder inside root with the same name.
    varieties (list of lists): varieties to consider for each region. The function will create a series of baselines for 
            each regions and each varieties for that region. The name of the variety will be used to name the baselines.
            Elements in the list of lists are string, for example [['arabica','robusta'],['arabica'],...]
    regionsIn (str): Name of the folder inside the regions folders where to find the input modis rasters to 
            be used for the averaging. These files are outputs of the smoothing. 
            It should be the same for all the regions 
    avgWeights (str): Raster (single layer) or shapefile with the weights to use for the averages. If shapefile, will be rasterized. 
            If raster, will be resampled to the correct resolution if needed. In each case the information in each 
            pixel/polygon should be a density of crop of interest (0 < d < 1).
    weightField (str): Only used if avgWeights is a shapefile. Name of the field containing the weight/density information
    startAvg (str): Starting date for the files to consider in the averaging (included). There will be one value per date.
    endAvg (str): Ending date for the files to consider in the averaging (included).
    alltouch (boolean): true or false, whether all the pixels touched by the region should be considered in the average or 
            only the pixels with their centroid inside the region.
            This parameter is used when rasterizing the shapefile with the weights and will therefore only be used if 
            avgWeights is a shapefile.
        
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
    
    # Create an empty dictionary to get the values for each of the regions
    averages = {}
    colnames = []
    
    for r in range(len(regions)):
        for v in range(len(varieties[r])):
            if not avgWeights[r][v]:
                continue
            
            colnames.append(regions[r] + '_' + varieties[r][v])
            
            # Get the images to consider
            print('Averaging region ' + str(regions[r]) + '...')
            
            # Import all the smooth modis images on disk
            onDisk = [f for f in os.listdir(os.path.join(root, regions[r], regionsIn)) if f.endswith('.tif')]
            
            # Dates of these files
            datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
            # Transform into date format
            datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
            # Keep only the files and dates within the dates to process
            onDisk = [f for f, d in zip(onDisk, datesAll) if d >= startAvg and d <= endAvg]
            datesAll = [d for d in datesAll if d >= startAvg and d <= endAvg]
            
            if not onDisk:
                print('no modis images to process in ' + regions[r])
                continue
            
            if avgWeights[r] and avgWeights[r][v]:
                avgW = avgWeights[r][v]
            else:
                avgW = None
            
            # Compute the averages fo rall the dates
            avgRegion = avgRegionRaster(images=onDisk,
                                        datesImg=datesAll,
                                        weightsRaster=avgW,
                                        weightField=weightField,
                                        alltouch=False,
                                        BlockXSize=256, BlockYSize=256)
            
            if not avgRegion:
                break
            
            # Transform the results into a dictionary easier to export
            for k, s in avgRegion.iteritems():
                # Transform into 'per Hectare'
                s = s / (250.) * 10000.
                if k in averages:
                    averages[k][regions[r] + '_' + varieties[r][v]] = s
                else:
                    averages[k] = {}
                    averages[k]['date'] = k
                    averages[k][regions[r] + '_' + varieties[r][v]] = s
            
    # Export the dictionary
    outNm = 'Weighted_avg_ndvi_' + startAvg.strftime('%Y-%m-%d') + '_' + endAvg.strftime('%Y-%m-%d') + '.txt'
    # Sort the dates
    datesAll.sort()
    # order the output by date in a list. Each element is an element of the original dictionary and will be exported
    out = []
    for date in datesAll:
        out.append(averages[date.strftime('%Y-%m-%d')])
    with open(os.path.join(root, outNm), "w") as f:
        dict_writer = DictWriter(f, ['date'] + colnames, extrasaction='ignore', delimiter="\t", restval="0")
        dict_writer.writeheader()
        for p in out:
            dict_writer.writerow(p)


def avgRegionRaster(images, datesImg, weightsRaster=None, weightField=None,
                    alltouch=False, BlockXSize=256, BlockYSize=256):
    
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
        weightsRaster = new_raster_from_base(baseImg, 'temp', 'MEM', -1, gdal.GDT_Float32) 
        
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
            weightsRaster = warp_raster(weightsRaster, baseImg, resampleOption='nearest', outputURI=None, outFormat='MEM')
    
    baseImg = None
    nodataWeights = weightsRaster.GetRasterBand(1).GetNoDataValue()
    
    # Prepare an empty dictionary to holde the results
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
        xBlocks = int(round(xsize / BlockXSize)) + 1
        yBlocks = int(round(ysize / BlockYSize)) + 1
        band = None
        
        sumNdvi = []
        sumWeights = []
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                # Read the block from the image
                blockBase = readRasterBlock(baseImg, xStep * BlockXSize, yStep * BlockYSize, BlockXSize, BlockYSize)
                
                # Read the block from the weights
                blockWeight = readRasterBlock(weightsRaster, xStep * BlockXSize, yStep * BlockYSize, BlockXSize, BlockYSize)
                
                # Recast the type to be sure
                blockBase = blockBase.astype(np.float32)
                blockWeight = blockWeight.astype(np.float32)
                
                # Replace the no data values by 0
                blockBase[np.logical_or(np.logical_or(blockBase == nodataBase, np.isnan(blockBase)), np.logical_or(blockBase > 1, blockBase < -1))] = 0.
                blockWeight[np.logical_or(np.logical_or(blockWeight == nodataWeights, np.isnan(blockWeight)), np.logical_or(blockWeight > 1, blockWeight < 0))] = 0.
            
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


def maskQualityVI(ndviRaster, qualityRaster, outRaster=None, nodata=None):
    '''
    Function to mask a composite modis ndvi image with a quality layer.
    
    Each pixel of the quality layer is converted to 16 bit to screen
    quality problems and mask pixels with low quality.
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
    if not nodata:
        nodata = ndvi.GetRasterBand(1).GetNoDataValue()
        if not nodata:
            nodata = 65535
                    
    # Loop through the arrays to mask the ndvi
    
    for p, q in np.nditer([ndviArray, qArray],
                          flags=['reduce_ok'],
                          op_flags=[['readwrite'], ['readonly']],
                          op_dtypes=[ndviType, qType]):
        p[...] = screenQualityVI(p, q, nodata)
    
    if outRaster:
        # Create an empty raster copy
        out = new_raster_from_base(ndvi, outRaster, 'GTiff', nodata,
                                   getattr(gdal,
                                           getGDALTypeFromNumber(
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
    
    
def screenQualityVI(pixel, index, nodata):
    # Function to screen the quality of a pixel of ndvi in composite 
    # modis product based on quality layer.
    
    if pixel == nodata or np.isnan(pixel):
        return nodata
    
    bit16 = np.binary_repr(index, width=16)
    
    if (bit16[2:5] == '001' and bit16[14:16] in ['00', '01', '10'] and 
        bit16[9:13] in ['0000', '0001', '0010', '0011', '0100', '0101', '0111',
                        '1000']):
         p = pixel
    else:
        p = nodata
    
    return(p)


def percMissingStack(images, outName, nodata=None):
    '''
    Takes a list of images/rasters of different dates and returns a 
    single raster wiwth for each pixel the % of missing
    
    A no data value should be specified if it is missing from the rasters
    The no data value should be the same for all the rasters
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
    out = new_raster_from_base(toProcess[0], outName, 'GTiff', np.nan, gdal.GDT_Float32)
    
    # Get the size of the rasters to identify the limits of the blocks to loop through
    band = toProcess[0].GetRasterBand(1)
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    # Set the block size
    BlockXSize = 256 * 4
    BlockYSize = 256 * 4
    # Get the number of blocks in x and y directions
    xBlocks = int(round(xsize / BlockXSize)) + 1
    yBlocks = int(round(ysize / BlockYSize)) + 1
    
    totBlocks = xBlocks * yBlocks
    progress = 1
        
    for xStep in range(xBlocks):
        for yStep in range(yBlocks):
            print('Processing block ' + str(progress) + ' of ' + str(totBlocks))
            progress += 1
            
            blocks = [readRasterBlock(p, xStep * BlockXSize, yStep * BlockYSize, BlockXSize, BlockYSize) for p in toProcess]
            
            # Bring the blocks together into one single array
            blocks = np.dstack(blocks)
            
            # Do the smoothing for each pixel
            blocks = percMissing(blocks, nodataImg)
            
            # Change the values in the output raster
            out.GetRasterBand(1).WriteArray(blocks, xStep * BlockXSize, yStep * BlockYSize)
        
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


def warp_raster(src, dst, resampleOption='nearest', outputURI=None, outFormat='MEM'):
    """
    ---------------------------------------------------------------------------------------------
    Function : Warp a source raster to the resolution, extent and projection of a destination raster.
           
            The function returns the resulting raster. If outFormat is different from 'MEM', the 
            raster is also saved to disk using the information provided in outputURI.
            
            Inputs
            --src (gdal Dataset): source raster to be warped
            --dst (gdal Dataset): destination raster that will provide the resolution, extent and projection
            --resampleOption (string): One of 'nearest', 'bilinear', 'cubic', 'cubic spline', 'lanczos', 
                    'average', or 'mode'. Method to use to resample the pixels of the source raster
            --outputURI (string, optional): Full address and name of the output raster. If outFormat is 'MEM',
                    this argument is ignored and the function simply produces a raster in memory. 
                    The extension for the output file should match the outFormat.
            
            --outFormat (string, optional): Format to use for the output raster from the function. 
                    Use 'GTiff' for a .tif output. Default creates a raster in memory.
    ---------------------------------------------------------------------------------------------
    """  
    if not type(src) is gdal.Dataset:
        return False
    
    if not type(dst) is gdal.Dataset:
        return False
    
    # Define resampling options
    resampleOptions = {'nearest': gdalconst.GRA_NearestNeighbour, 'bilinear':gdalconst.GRA_Bilinear,
                   'cubic':gdalconst.GRA_Cubic, 'cubic spline':gdalconst.GRA_CubicSpline,
                   'lanczos':gdalconst.GRA_Lanczos, 'average':gdalconst.GRA_Average,
                   'mode':gdalconst.GRA_Mode} 
    
    if not resampleOption in resampleOptions.keys():
        return False
    
    # Raster to host the warped output
    if outFormat == 'MEM':
        rOut = new_raster_from_base(dst, 'temp', 'MEM', 0, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    else:
        rOut = new_raster_from_base(dst, outputURI, outFormat, 0, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    
    # Warp: the parameters are source raster, destination raster, source projection, destination projection, resampling option 
    gdal.ReprojectImage(src, rOut, src.GetProjection(), rOut.GetProjection(), resampleOptions[resampleOption])
    
    return rOut


def readRasterBlock(src, xStart, yStart, xBlockSize, yBlockSize, band=1):
    '''
    Function to read a block of data from a gdal Dataset. It will return the block as a numpy array
    
    src (gdal Dataset): raster to extract the block from
    xStart (int): X pixel value from which to start extraction (upper left corner)
    yStart (int): Y pixel value from which to start extraction (upper left corner)
    xBlockSize (int): X size of the block to extract. If the block is larger than the number 
            of pixels from xStart, will extract up to the raster limit
    yBlockSize (int): Y size of the block to extract. If the block is larger than the number 
            of pixels from xStart, will extract up to the raster limit
    band (int): band from which to extract the block. 1 by default.
    '''
    
    # Get the wanted band
    band = src.GetRasterBand(band)
    
    # Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    
    if yStart + yBlockSize < ysize:
        rows = yBlockSize
    else:
        rows = ysize - yStart
    
    if xStart + xBlockSize < xsize:
        cols = xBlockSize
    else:
        cols = xsize - xStart
            
    outArray = band.ReadAsArray(xStart, yStart, cols, rows)
    
    return outArray


def new_raster_from_base(base, outputURI, formatR, nodata, datatype, bands=None):
    ''' 
    Create an empty copy of a raster from an existing one
    
    base: gdal raster layer
        Name of the variable with the input raster to copy
    
    outputURI: string
        Address + name of the output raster (extension should agree with format, none for memory)
        
    formatR: string
        Format for the dataset (e.g. "GTiff", "MEM")
    
    nodata: int/float
        No data value (type should agree with raster type)
    
    datatype: gdal data type (e.g. gdal.GDT_Int32)
        Data type for the raster
    
    bands: [optional] int
        Number of bands for the output raster. 
        If not specified, will use the number of bands of the input raster
    
    The function returns a gdal raster variable filled with the nodata value
    '''
    
    cols = base.RasterXSize
    rows = base.RasterYSize
    projection = base.GetProjection()
    geotransform = base.GetGeoTransform()
    if not bands:
        bands = base.RasterCount

    driver = gdal.GetDriverByName(formatR)
    
    if formatR == "GTiff":
        new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype, options=['COMPRESS=LZW', 'BIGTIFF=IF_NEEDED'])
    else:
        new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype)
    new_raster.SetProjection(projection)
    new_raster.SetGeoTransform(geotransform)

    for i in range(bands):
        new_raster.GetRasterBand(i + 1).SetNoDataValue(nodata)
        new_raster.GetRasterBand(i + 1).Fill(nodata)

    return new_raster


def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    Returns values based on origin and resolution of the raster. 
    The actual values can fall outside of the raster and should be checked.
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    col = int((x - ulX) / xDist)
    row = int((ulY - y) / xDist)
    
    return (col, row)


def getGDALTypeFromNumber(nb):
    if nb == 0:
        out = 'GDT_Unknown'
    elif nb == 1:
        out = 'GDT_Byte'
    elif nb == 2:
        out = 'GDT_UInt16'
    elif nb == 3:
        out = 'GDT_Int16'
    elif nb == 4:
        out = 'GDT_UInt32'
    elif nb == 5:
        out = 'GDT_Int32'
    elif nb == 6:
        out = 'GDT_Float32'
    elif nb == 7:
        out = 'GDT_Float64'
    elif nb == 8:
        out = 'GDT_CInt16'
    elif nb == 9:
        out = 'GDT_CInt32'
    elif nb == 10:
        out = 'GDT_CFloat32'
    elif nb == 11:
        out = 'GDT_CFloat64'
    elif nb == 12:
        out = 'GDT_TypeCount'
    
    return(out)
