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
    7- Compute reference rasters with decile values
    8- Compute rasters showing position of each pixel compared to reference deciles
    9- Compute averages by AOI
    
"""

import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

import os, re, multiprocessing
import pathos.multiprocessing as pmp
import MODIS_gedata_toolbox as md  # GEDATA toolbox for MODIS related tools
import gapfill  # Python implementation of the interpolation algorithm
from datetime import datetime
from csv import DictWriter, writer


def main():
    #####################################################################################################################
    #####################################################################################################################
    # #PARAMETERS
    
    # Allow parallel computing?
    allowPara = True
    # Number of cores to use?
    nCores = 3
    
    # Satellite
    satelliteModis = 'terra'  # 'terra' # 'aqua'
    
    # Root folder
    prefixRootSys = '/media/olivier/olivier_ext1/gedata_current/jde_coffee'  #'E:/gedata_current' #'/home/olivierp/jde_coffee'
    
    # #DIRECTORIES parameters
    # Working directory
    dst = os.path.join(prefixRootSys, 'MODIS/collection6/' + satelliteModis + '/Brazil')
    # Directory data sources
    dataDir = os.path.join(prefixRootSys, 'data/Brazil')
    # Folder inside dst to use for temporary files (should be empty)
    tempDir = os.path.join(prefixRootSys, 'Temp')
    # Destination folder for the downloaded MODIS tiles
    rawdataDir = os.path.join(dst, 'raw_data')
    
    # #REGIONS parameters
    # Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    # Names of the regions (also folders names) 
    states = ["CER", "CHA", "CO", "ES", "MO", "SDM", "SP", "ZM"]
    # Varieties in each case
    varieties = [['arabica'], ['arabica'], ['arabica'], ['arabica', 'robusta'],
                 ['arabica'], ['arabica'], ['arabica'], ['arabica', 'robusta']]
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
    fillMissing = True
    # Smooth images
    smooth = False
    # Create baselines
    createBaselines = False
    # Rank individual images against baseline images
    ranking = False
    # Average MODIS value in each region
    avgValue = False
    #Additional options with the average
    # Compute the quality information to inform the 
    # average values for each region
    qualIndex = True
    #Create charts of the differential to the long term
    chartDiff = True
    
    ########### DOWNLOAD
    # Product to download
    if satelliteModis == 'terra':
        product = 'MOD13Q1.006'
    else:
        product = 'MYD13Q1.006'
    # Username for the earthdata website
    user = "olivierpfrancois"
    # Password for the earthdata website
    pwd = "Michele1950"
    # Tiles to download.
    tiles = ['h13v10', 'h13v11', 'h14v10', 'h14v11']  # ['h28v07']
    # Start date for the product download (format YYYY-MM-DD)
    #    If None, will default to date of most recent MODIS file on disk if any, or stop the process
    startDownload = '2018-03-01'
    # End date for the product download (format YYYY-MM-DD)
    #    If None, defaults to today
    endDownload = None
    
    ############ MOSAIC
    # Starting date for the files to mosaic
    #    If None, will default to the files that have been just downloaded if 
    #    any.
    startMosaic = '2018-03-01'  # '2017-02-01'
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
    startCheck = '2018-03-01'
    # End date for the files to check
    endCheck = None
    
    ############ FILL MISSING
    # Input folder for the images to fill
    inMissing = statesMaskedFolder
    # Output folder for the images to fill
    outMissing = statesFilledFolder
    # Year(s) of images to fill
    yearsMissing = [2017,2018]
    # Day(s) of images to fill
    daysMissing = None #[[49,65,81]]
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
    startSmooth = '2016-06-15'  # '2012-03-01'
    # Ending date for the files to include as input in the smoothing process
    #    If None, defaults to today
    endSmooth = None
    # Start and end dates for the files to save to the disk after smoothing
    #    If None, defaults to 6 months before end smoothing date
    startSaveS = '2017-08-01'
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
    maskBaseline = [[dst + '/CER/' + 'masks/CER_densities_arabica_from_classifications_250m.tif'],
                    [dst + '/CHA/' + 'masks/CHA_densities_arabica_from_classifications_250m.tif'],
                    [dst + '/CO/' + 'masks/CO_densities_arabica_from_classifications_250m.tif'],
                    [dst + '/ES/' + 'masks/ES_densities_arabica_from_classifications_250m.tif',
                            dst + '/ES/' + 'masks/ES_densities_robusta_from_classifications_250m.tif'],
                    [dst + '/MO/' + 'masks/MO_densities_arabica_from_classifications_250m.tif'],
                    [dst + '/SDM/' + 'masks/SDM_densities_arabica_from_classifications_250m.tif'],
                    [dst + '/SP/' + 'masks/SP_densities_arabica_from_classifications_250m.tif'],
                    [dst + '/ZM/' + 'masks/ZM_densities_arabica_from_classifications_250m.tif',
                            dst + '/ZM/' + 'masks/ZM_densities_robusta_from_classifications_250m.tif']]
    outModelRaster = [[dst + '/CER/' + 'masks/CER_densities_arabica_from_classifications_1km.tif'],
                      [dst + '/CHA/' + 'masks/CHA_densities_arabica_from_classifications_1km.tif'],
                      [dst + '/CO/' + 'masks/CO_densities_arabica_from_classifications_1km.tif'],
                      [dst + '/ES/' + 'masks/ES_densities_arabica_from_classifications_1km.tif',
                            dst + '/ES/' + 'masks/ES_densities_robusta_from_classifications_1km.tif'],
                      [dst + '/MO/' + 'masks/MO_densities_arabica_from_classifications_1km.tif'],
                      [dst + '/SDM/' + 'masks/SDM_densities_arabica_from_classifications_1km.tif'],
                      [dst + '/SP/' + 'masks/SP_densities_arabica_from_classifications_1km.tif'],
                      [dst + '/ZM/' + 'masks/ZM_densities_arabica_from_classifications_1km.tif',
                            dst + '/ZM/' + 'masks/ZM_densities_robusta_from_classifications_1km.tif']]
    
    ############ RANKING 
    # Ranking individual dates modis images in terms of deciles using 
    #    the baselines
    # Starting and ending dates for the images to consider. Included
    #   If None, will default to 60 days before the endRank date
    startRank = '2017-12-15'
    #   If None, will default to today
    endRank = None
    # Minimum density of coffee to consider 
    minCoffee = [0.05, 0.15]
    # File to use for masking the output using the density
    # maskRank = [['masks/LD_densities_coffee_from_classifications_1km.tif']]
    maskRank = [[dst + '/CER/' + 'masks/CER_densities_arabica_from_classifications_1km.tif'],
                [dst + '/CHA/' + 'masks/CHA_densities_arabica_from_classifications_1km.tif'],
                [dst + '/CO/' + 'masks/CO_densities_arabica_from_classifications_1km.tif'],
                [dst + '/ES/' + 'masks/ES_densities_arabica_from_classifications_1km.tif',
                    dst + '/ES/' + 'masks/ES_densities_robusta_from_classifications_1km.tif'],
                [dst + '/MO/' + 'masks/MO_densities_arabica_from_classifications_1km.tif'],
                [dst + '/SDM/' + 'masks/SDM_densities_arabica_from_classifications_1km.tif'],
                [dst + '/SP/' + 'masks/SP_densities_arabica_from_classifications_1km.tif'],
                [dst + '/ZM/' + 'masks/ZM_densities_arabica_from_classifications_1km.tif',
                    dst + '/ZM/' + 'masks/ZM_densities_robusta_from_classifications_1km.tif']]
    
    ############ AVERAGE
    # Starting and ending dates for the images to consider. Included
    #    If None, defaults to 1 year before today
    startAvg = '2006-01-01'
    #    If None, will default to today
    endAvg = None
    #Grid/Raster to use for the averaging --> FULL PATH!!!!!!!
    avgWeights = [[dst + '/CER/masks/CER_densities_arabica_from_classifications_250m.tif'],
                  [dst + '/CHA/masks/CHA_densities_arabica_from_classifications_250m.tif'],
                  [dst + '/CO/masks/CO_densities_arabica_from_classifications_250m.tif'],
                  [dst + '/ES/masks/ES_densities_arabica_from_classifications_250m.tif',
                   dst + '/ES/masks/ES_densities_robusta_from_classifications_250m.tif'],
                  [dst + '/MO/masks/MO_densities_arabica_from_classifications_250m.tif'],
                  [dst + '/SDM/masks/SDM_densities_arabica_from_classifications_250m.tif'],
                  [dst + '/SP/masks/SP_densities_arabica_from_classifications_250m.tif'],
                  [dst + '/ZM/masks/ZM_densities_arabica_from_classifications_250m.tif',
                   dst + '/ZM/masks/ZM_densities_robusta_from_classifications_250m.tif']]
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
                               pwd=pwd, user=user, tiles=tiles, product=product,
                               startDownload=startDownload,
                               endDownload=endDownload,
                               satellite=satelliteModis)
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
        print('Starting quality check')
        
        if not endCheck:
            endCheck = datetime.now()
            endCheck = endCheck.date()
        else:
            endCheck = datetime.strptime(endCheck, '%Y-%m-%d').date()
        
        if startCheck:
            startCheck = datetime.strptime(startCheck, '%Y-%m-%d').date()
        
        for s, b in zip(states, statesBoundFiles):
            print('   Starting region ' + s)
            
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
                p = pmp.Pool(nCores)
                
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
                
                p.close()
                p.join()
                
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
        print('Starting interpolation of missing values')
        
        for s, b in zip(range(len(states)), statesBoundFiles):
            print('   Starting region ' + states[s])
            if not s>4:
                continue
            inputRasters = [os.path.join(dst, states[s], inMissing, f) for 
                            f in os.listdir(os.path.join(dst, states[s], inMissing)) 
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
            '''
            expans = gapfill.gapFill(rasters=inputRasters, seasons=days, years=years,
                                     outFolder=os.path.join(dst, states[s], outMissing),
                                     suffix=suffMissing, nodata=[-3000], iMax=20,
                                     subsetSeasons=daysMissing, subsetYears=yearsMissing, 
                                     subsetMissing=None, clipRange=(-2000, 10000), 
                                     parallel=allowPara, nCores=nCores)
            '''
            if avgWeights[s]:
                avgW = avgWeights[s]
            else:
                avgW = None
            expans = gapfill.gapFillTest(rasters=inputRasters, seasons=days, years=years,
                            outFolder=os.path.join(dst, states[s], outMissing),
                            suffix=suffMissing, nodata=[-3000], iMax=27,
                            subsetSeasons=daysMissing, subsetYears=yearsMissing, subsetMissing=None,
                            clipRange=(-2000, 10000), maskRaster=avgW, 
                            parallel=allowPara, nCores=nCores)
            
            with open(dst+"/test"+states[s]+".csv", "wb") as f:
                exp = writer(f)
                exp.writerows(expans)
                
            # Mask the resulting rasters to the specific extent of the AOI
            for r, y, d in zip(inputRasters, years, days):
                if (yearsMissing and not y in yearsMissing):
                    continue
                elif (yearsMissing and daysMissing and 
                      not d in daysMissing[yearsMissing.index(y)]):
                        continue
                elif not yearsMissing and daysMissing and not d in daysMissing:
                    continue
                
                nameR = os.path.join(dst, states[s], outMissing,
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
        
        # Create an empty dictionary to get the values for each of the regions
        averages = {}
        if qualIndex:
            qualityMasked = {}
            qualityFilled = {}
        colnames = []
        datesAll = []
        
        for r in range(len(states)):
            for v in range(len(varieties[r])):
                if not avgWeights[r][v]:
                    continue
                
                colnames.append(states[r] + '_' + varieties[r][v])
                
                # Get the images to consider
                print('Averaging region ' + str(states[r]) + '...')
                
                if avgWeights[r] and avgWeights[r][v]:
                    avgW = avgWeights[r][v]
                else:
                    avgW = None
                
                avgRegion = md.avgRegionRasterWrap(
                    regionIn=os.path.join(dst, states[r], statesSmoothFolder),
                    avgWeights=avgW,
                    weightField=weightField,
                    startAvg=startAvg, endAvg=endAvg,
                    alltouch=False)
                '''
                avgRegion = md.avgRegionQualWrap(
                    regionIn=os.path.join(dst, states[r], statesSmoothFolder),
                    maskedIn=os.path.join(dst, states[r], statesMaskedFolder),
                    avgWeights=avgW,
                    weightField=weightField,
                    startAvg=startAvg, endAvg=endAvg,
                    alltouch=False)
                '''
                if avgRegion == False:
                    break
                if len(avgRegion) == 0:
                    continue
                
                if qualIndex:
                    # Compute the quality for all the dates
                    qualMasked = md.computeQualityIndexNdviWrap(
                        regionIn=os.path.join(dst, states[r], statesMaskedFolder),
                        avgWeights=avgW,
                        weightField=weightField,
                        missingValue=-3000,
                        startAvg=startAvg, endAvg=endAvg,
                        alltouch=False)
                    
                    qualFilled = md.computeQualityIndexNdviWrap(
                        regionIn=os.path.join(dst, states[r], statesFilledFolder),
                        avgWeights=avgW,
                        weightField=weightField,
                        missingValue=None,
                        startAvg=startAvg, endAvg=endAvg,
                        alltouch=False)
                    
                # Transform the results into a dictionary easier to export
                for k, s in avgRegion.iteritems():
                    # Transform into 'per Hectare'
                    s = s / (250.) * 10000.
                    if not k in datesAll:
                        datesAll.append(k)
                    if not k in averages:
                        averages[k] = {}
                        averages[k]['date'] = k
                    averages[k][states[r] + '_' + varieties[r][v]] = s
                        
                    if qualIndex:
                        if qualMasked:
                            if not k in qualityMasked:
                                qualityMasked[k] = {}
                                qualityMasked[k]['date'] = k
                            if k in qualMasked.keys():
                                qualityMasked[k][states[r] + '_' + varieties[r][v]] = qualMasked[k]
                            else:
                                qualityMasked[k][states[r] + '_' + varieties[r][v]] = ' '
                        
                        if qualFilled:
                            if not k in qualityFilled:
                                qualityFilled[k] = {}
                                qualityFilled[k]['date'] = k
                            if k in qualFilled.keys():
                                qualityFilled[k][states[r] + '_' + varieties[r][v]] = qualFilled[k]
                            else:
                                qualityFilled[k][states[r] + '_' + varieties[r][v]] = ' '
        
        # Sort the dates and get min and max
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        datesAll.sort()
        outMin = min(datesAll)
        outMax = max(datesAll)
        
        # Create output name
        outNm = ('Weighted_avg_ndvi_' + outMin.strftime('%Y-%m-%d') + '_' + 
                 outMax.strftime('%Y-%m-%d') + '.txt')
        # order the output by date in a list. Each element is an element of the original dictionary and will be exported
        out = []
        for date in datesAll:
            out.append(averages[date.strftime('%Y-%m-%d')])
        # Export the dictionary
        with open(os.path.join(dst, outNm), "w") as f:
            dict_writer = DictWriter(f, ['date'] + colnames, extrasaction='ignore', delimiter="\t", restval="0")
            dict_writer.writeheader()
            for p in out:
                dict_writer.writerow(p)
        
        # Same process for the quality information
        if qualIndex:
            if qualityMasked:
                outNmMasked = ('Quality_masked_ndvi_' + outMin.strftime('%Y-%m-%d') + '_' + 
                               outMax.strftime('%Y-%m-%d') + '.txt')
                outM = []
                for date in datesAll:
                    outM.append(qualityMasked[date.strftime('%Y-%m-%d')])
                with open(os.path.join(dst, outNmMasked), "w") as f:
                    dict_writer = DictWriter(f, ['date'] + colnames, extrasaction='ignore', delimiter="\t", restval="0")
                    dict_writer.writeheader()
                    for p in outM:
                        dict_writer.writerow(p)
                
            if qualityFilled:
                outNmFilled = ('Quality_filled_ndvi_' + outMin.strftime('%Y-%m-%d') + '_' + 
                               outMax.strftime('%Y-%m-%d') + '.txt')
                outF = []
                for date in datesAll:
                    outF.append(qualityFilled[date.strftime('%Y-%m-%d')])
                with open(os.path.join(dst, outNmFilled), "w") as f:
                    dict_writer = DictWriter(f, ['date'] + colnames, extrasaction='ignore', delimiter="\t", restval="0")
                    dict_writer.writeheader()
                    for p in outF:
                        dict_writer.writerow(p)

        if chartDiff:
            #Remove the date from the variables in averages
            for v in averages.itervalues():
                v.pop('date', None)
            #Create the plots
            md.plotModisLtavg(inDic=averages, ltAvgStart=2006, ltAvgEnd=2016, 
                              dateStartChart='07-01', yearsPlot=range(2009,2019), 
                              outFolder=dst)
            
if __name__ == '__main__':
    main()
