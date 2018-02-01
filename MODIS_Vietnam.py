# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Script to download MODIS images from the server, 
    and mosaic the images and mask them for each region
    to process.
    
    The image treatment process for each AOI is the following:
    1- Download the MODIS tiles
    2- Mosaic and cut the the extent of each AOI
    3- Mask missing/bad quality pixels in entire extent and mark those that 
        need to be filled
    4- Fill missing values
    5- Mask the image to only keep the pixels of interest
    6- Smooth the images temporally
    
"""

import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

import os, re, multiprocessing
import pathos.multiprocessing as mp
import MODIS_gedata_toolbox as md  # GEDATA toolbox for MODIS related tools
import gapfill  # Python implementation of the interpolation algorithm
from datetime import datetime, timedelta


def main():
    #####################################################################################################################
    #####################################################################################################################
    # #PARAMETERS
    
    # Allow parallel computing?
    allowPara = True
    # Number of cores to use?
    nCores = 3
    
    # Root folder
    prefixRootSys = '/media/olivier/olivier_ext1/gedata_current/jde_coffee'  # 'E:/gedata_current' # '/home/olivierp' 
    
    # #DIRECTORIES parameters
    # Working directory
    dst = os.path.join(prefixRootSys, 'MODIS/collection6/Vietnam')
    # Directory data sources
    dataDir = os.path.join(prefixRootSys, 'data/Vietnam')
    # Folder inside dst to use for temporary files (should be empty)
    tempDir = os.path.join(prefixRootSys, 'Temp')
    # Destination folder for the downloaded MODIS tiles
    rawdataDir = os.path.join(dst, 'raw_data')
    
    # #REGIONS parameters
    # Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    # Names of the regions (also folders names) 
    states = ["LD"]
    # Varieties in each case
    varieties = [['coffee']]
    # Addresses of the shapefiles with the boundaries for each of the regions
    # !!!!SHOULD BE IN EPSG 4326 PROJECTION
    statesBoundFiles = [dataDir + '/' + s + '/aoi/AOI_' + s + '.shp' for s in states] 
    # Name of subfolder in each region where to save the raw mosaic data 
    # (should be in the folders of the regions)
    statesRawFolder = 'raw_data'
    # Name of the subfolder where to save the masked images
    statesMaskedFolder = 'masked_missing'
    # Name of the subfolder where to save the filled images
    statesFilledFolder = 'filled_missing'
    # Name of subfolder where to save the smoothed mosaic data (should be in the folders of the regions)
    statesSmoothFolder = 'smooth_data'
    # Name of subfolder where to save (if produced) and get the reference 
    # baseline for each date in the year (should be in the folders of the regions)
    statesRefFolder = 'baseline'
    
    # Process decision dummies
    
    # Download images
    dload = False
    # Mosaic images for each region and crop to extent 
    mosaic = False
    # Check quality 
    checkQuality = False
    # Fill missing values and mask by exact AOI 
    fillMissing = False
    # Smooth images
    smooth = False
    # Create baselines
    createBaselines = False
    # Rank individual images against baseline images
    ranking = False
    # Average MODIS value in each region
    avgValue = False
    
    # Get the percentage of missing values for each pixel across history
    checkPercMissing = False
    maskLand = False
    
    ########### DOWNLOAD
    # Product to download
    product = 'MOD13Q1.006'
    # Username for the earthdata website
    user = "olivierpfrancois"
    # Password for the earthdata website
    pwd = "Michele1950"
    # Tiles to download.
    tiles = ['h28v07']  # ['h28v07']
    # Start date for the product download (format YYYY-MM-DD)
    #    If None, will default to date of most recent MODIS file on disk if any, or stop the process
    startDownload = '2016-01-01'  # '2017-05-26'
    # End date for the product download (format YYYY-MM-DD)
    #    If None, defaults to today
    endDownload = None
    
    ############ MOSAIC
    # Starting date for the files to mosaic
    #    If None, will default to the files that have been just downloaded if 
    #    any.
    startMosaic = '2016-01-01'
    # startMosaic = '2005-01-01'
    # Ending date for the files to mosaic
    #    If None, defaults to today
    endMosaic = None  # '2005-02-01'
    
    ############ MASK QUALITY
    # Input folder of the images to mask  
    inCheck = statesRawFolder
    # Output folder of the images to mask
    outCheck = statesMaskedFolder
    # Start date for the files to check
    startCheck = '2016-01-01'
    # End date for the files to check
    endCheck = None
    
    ############ FILL MISSING
    # Input folder for the images to fill
    inMissing = statesMaskedFolder
    # Output folder for the images to fill
    outMissing = statesFilledFolder
    # Year(s) of images to fill
    yearsMissing = [2016, 2017, 2018]
    # Day(s) of images to fill
    daysMissing = None
    # Suffix to put at the end of the name of the 
    # images after filling
    suffMissing = 'f'
    # !!! The two conditions are additive (AND)
    
    ############ SMOOTH
    # Size of the regression and average window (Swets algorithm)
    regWindow = 7
    avgWindow = 3
    # Starting date for the files to include as input in the smoothing process
    #    If None, defaults to 1 year before the end smoothing date
    startSmooth = '2015-01-01'  # '2012-03-01'
    # Ending date for the files to include as input in the smoothing process
    #    If None, defaults to today
    endSmooth = None
    # Start and end dates for the files to save to the disk after smoothing
    #    If None, defaults to 6 months before end smoothing date
    startSaveS = '2016-01-01'
    # startSaveS = '2017-03-01'
    endSaveS = None  # None to save them up to the end smoothing date
    
    ############ BASELINE
    # Starting and ending years for the reference period. Included.
    startRef = 2006
    endRef = 2016
    # Mask and output model to use for the baseline rasters. 
    # The mask should overlap perfectly with the modis images
    # The output model can have a different resolution, in which case the baseline will be produced with that resolution
    # maskBaseline = [['masks/LD_densities_coffee_from_classifications_250m.tif']]
    # outModelRaster = [['masks/LD_densities_coffee_from_classifications_1km.tif']]
    maskBaseline = [[dst + '/LD/' + 'masks/LD_densities_coffee_from_classifications_250m.tif']]
    outModelRaster = [[dst + '/LD/' + 'masks/LD_densities_coffee_from_classifications_1km.tif']]
    
    ############ RANKING 
    # Ranking individual dates modis images in terms of deciles using 
    #    the baselines
    # Starting and ending dates for the images to consider. Included
    #   If None, will default to 60 days before the endRank date
    startRank = '2016-06-01'
    #   If None, will default to today
    endRank = None
    # Minimum density of coffee to consider 
    minCoffee = [0.05, 0.15]
    # File to use for masking the output using the density
    # maskRank = [['masks/LD_densities_coffee_from_classifications_1km.tif']]
    maskRank = [[dst + '/LD/' + 'masks/LD_densities_coffee_from_classifications_1km.tif']]
    
    ############ AVERAGE
    # Starting and ending dates for the images to consider. Included
    #    If None, defaults to 1 year before today
    startAvg = '2006-01-01'
    #    If None, will default to today
    endAvg = None
    #Grid/Raster to use for the averaging --> FULL PATH!!!!!!!
    avgWeights = [[dst + '/LD/masks/LD_densities_coffee_from_classifications_250m.tif']]
    # Name of the field with the density information if the masks for averaging 
    #    are shapefiles 
    weightField = None
    
    #####################################################################################################################
    #####################################################################################################################
    # #ACTIVE CODE
    
    if allowPara and not nCores:
        nCores = multiprocessing.cpu_count()
    
    if dload:
        newHdf = md.downloadMODIS(dstFolder=os.path.join(dst, rawdataDir),
                               pwd=pwd, user=user, tiles=tiles, product=product, startDownload=startDownload,
                               endDownload=endDownload)
    else:
        newHdf = []
    
    if mosaic:
        
        thereHdf = [f for f in os.listdir(os.path.join(dst, rawdataDir)) if f.endswith('.hdf')]
        
        if not newHdf and not thereHdf:
            return
        
        else:
            if not startMosaic:
                print('No start date provided for mosaic. Will mosaic ' + 
                      'downloaded images only')
            print('Starting the mosaic process')
            
            md.mosaicMODISWrapper(root=dst,
                                  srcFolder=os.path.join(dst, rawdataDir),
                                  tmpFolder=tempDir,
                                  regions=states,
                                  regionsOut=statesRawFolder,
                                  regionsBoundaries=statesBoundFiles,
                                  tiles=tiles,
                                  subset='1 0 1 0 0 0 0 0 0 0 0 0',
                                  suffix=['NDVI', 'Quality'],
                                  nodataOut=[32767, 65535],
                                  startMosaic=startMosaic,
                                  endMosaic=endMosaic)
    
    if checkQuality:
        if not endCheck:
            endCheck = datetime.now()
            endCheck = endCheck.date()
        else:
            endCheck = datetime.strptime(endCheck, '%Y-%m-%d').date()
        
        if startCheck:
            startCheck = datetime.strptime(startCheck, '%Y-%m-%d').date()
        
        for s, b in zip(states, statesBoundFiles):
            # Import all the raw ndvi modis images on disk
            allNDVI = [os.path.join(dst, s, inCheck, f) for 
                        f in os.listdir(os.path.join(dst, s, inCheck)) 
                        if f.endswith('NDVI.tif')]
        
            # Dates of these files
            datesNDVI = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                        for f in allNDVI]
            # Transform into date format
            datesNDVI = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesNDVI]
            
            if not startCheck:
                startCheck = min(datesNDVI)
            
            # Keep only the files and dates within the dates to process
            allNDVI = [f for f, d in zip(allNDVI, datesNDVI) 
                       if d >= startCheck and d <= endCheck]
            allNDVI.sort()
            
            # Update the dates in the final selection
            datesNDVI = [d for d in datesNDVI if d >= startCheck and d <= endCheck]
            datesNDVI.sort()
            
            # Import all the quality modis images on disk
            allQuality = [os.path.join(dst, s, inCheck, f) for 
                          f in os.listdir(os.path.join(dst, s, inCheck)) 
                          if f.endswith('Quality.tif')]
            
            # Dates of these files
            datesQual = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) 
                        for f in allQuality]
            # Transform into date format
            datesQual = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesQual]
            
            # Keep the same files as the NDVI
            allQuality = [f for f, d in zip(allQuality, datesQual) if
                         d in datesNDVI]
            allQuality.sort()
            
            # Prepare the output names
            allOut = [re.sub(inCheck, outCheck, f) for f in allNDVI]
            
            # Define the dataset
            dataset = zip(allNDVI, allQuality, allOut)
            
            if allowPara:
                p = mp.Pool(nCores)
                
                p.map(lambda d: md.maskQualityVI(ndviRaster=d[0],
                                                 qualityRaster=d[1],
                                                 outRaster=d[2],
                                                 nodataOut=-3000),
                      dataset)
                
                p.map(lambda d: md.clipMaskRasterByShp(shp=b,
                                    raster=d[2],
                                    outRaster=d[2],
                                    clipR=False,
                                    maskR=True,
                                    dataToMask=[-3000],
                                    nodataOut=32767),
                      dataset)
            
            else:
                for d in dataset:
                    # Mask the low quality pixels
                    md.maskQualityVI(ndviRaster=d[0],
                                     qualityRaster=d[1],
                                     outRaster=d[2],
                                     nodataOut=-3000)
                    
                    # Mask the pixels outside the shapefile
                    md.clipMaskRasterByShp(shp=b,
                                        raster=d[2],
                                        outRaster=d[2],
                                        clipR=False,
                                        maskR=True,
                                        dataToMask=[-3000],
                                        nodataOut=32767)
    
    if fillMissing:
        
        for s, b in zip(states, statesBoundFiles):
            inputRasters = [os.path.join(dst, s, inMissing, f) for 
                            f in os.listdir(os.path.join(dst, s, inMissing)) 
                            if f.endswith('.tif')]
            
            inputRasters.sort()
            
            # Get the dates
            datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in inputRasters]
            # Transform into date format
            datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
            # Transform into days from start of the year
            days = [int(d.strftime('%j')) for d in datesAll]
            
            # Get the years for the files on disk
            years = [int(d.strftime('%Y')) for d in datesAll] 
            
            gapfill.gapFill(rasters=inputRasters, seasons=days, years=years,
                            outFolder=os.path.join(dst, s, outMissing),
                            suffix=suffMissing, nodata=[-3000], iMax=20,
                            subsetSeasons=daysMissing, subsetYears=yearsMissing, subsetMissing=None,
                            clipRange=(-2000, 10000), parallel=allowPara, nCores=nCores)
            
            # Mask the resulting rasters to the specific extent of the AOI
            for r, y, d in zip(inputRasters, years, days):
                if (yearsMissing and not y in yearsMissing):
                    continue
                elif (yearsMissing and daysMissing and 
                      not d in daysMissing[yearsMissing.index(y)]):
                        continue
                elif not yearsMissing and daysMissing and not d in daysMissing:
                    continue
                
                nameR = os.path.join(dst, s, outMissing,
                                     re.sub('.tif', '_' + suffMissing + '.tif',
                                            os.path.basename(r)))
                
                md.clipMaskRasterByShp(shp=b,
                                    raster=nameR,
                                    outRaster=nameR,
                                    clipR=False,
                                    maskR=True,
                                    dataToMask=None,
                                    nodataOut=32767)
        
    '''
    for s, b in zip(states, statesBoundFiles):
        inputRasters = [os.path.join(dst, s, outMissing, f) for 
                            f in os.listdir(os.path.join(dst, s, outMissing)) 
                            if f.endswith('.tif') and not '_f' in f]
        
        for r in inputRasters:
            md.clipMaskRasterByShp(shp=b,
                                    raster=r,
                                    outRaster=r, 
                                    clipR=False, 
                                    maskR=True, 
                                    dataToMask=None, 
                                    nodataOut=32767)
    '''
    if smooth:
        print('Starting the smoothing process')
        
        md.smoothMODISWrapper(root=dst,
                       regions=states,
                       regionsIn=statesFilledFolder,
                       regionsOut=statesSmoothFolder,
                       startSmooth=startSmooth,
                       endSmooth=endSmooth,
                       regWindow=regWindow,
                       avgWindow=avgWindow,
                       startSaveSmooth=startSaveS,
                       endSaveSmooth=endSaveS,
                       parallel=allowPara,
                       nCores=nCores)
    
    if createBaselines:
        print('Starting creation of baselines')
        md.createBaseline(root=dst, regions=states, varieties=varieties,
                          regionsIn=statesSmoothFolder,
                          regionsOut=statesRefFolder,
                          startRef=startRef, endRef=endRef, mask=maskBaseline,
                          outModelRaster=outModelRaster,
                          parallel=allowPara,
                          nCores=nCores)
    
    if ranking:
        print('Starting analysis of individual dates compared to baseline')
        md.rankDatesDeciles(root=dst, regions=states, varieties=varieties,
                            regionsIn=statesSmoothFolder,
                            refDecilesIn=statesRefFolder,
                            startRank=startRank, endRank=endRank,
                            mask=maskRank, minDensity=minCoffee)

    if avgValue:
        print('Computing the average ndvi value for each zone')
        md.computeAvgNdvi(root=dst, regions=states, varieties=varieties,
                          regionsIn=statesSmoothFolder, avgWeights=avgWeights,
                          weightField=weightField, startAvg=startAvg,
                          endAvg=endAvg)

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


if __name__ == '__main__':
    main()

