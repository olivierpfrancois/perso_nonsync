# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Script to download MODIS images from the server, 
    and mosaic the images and mask them for each region
    to process.
 """
import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

import os, re
import MODIS_gedata_toolbox as md
import gapfill
from datetime import datetime

def main():
    #####################################################################################################################
    #####################################################################################################################
    # #PARAMETERS
    
    prefixRootSys = '/media/olivier/olivier_ext1/gedata_current' # 'E:/gedata_current' # '/home/olivierp' 
    
    # #DIRECTORIES parameters
    # Working directory
    dst = prefixRootSys + '/jde_coffee/MODIS/collection6/Brazil'
    # Directory data sources
    origin = prefixRootSys + '/jde_coffee/data/Brazil'
    # Folder inside dst to use for temporary files (should be empty)
    tempDir = prefixRootSys + '/jde_coffee/Temp'
    # 'E:/gedata_current/jde_coffee/Temp'
    # Destination folder for the download
    rawdataDir = 'raw_data'
    
    # #REGIONS parameters
    # Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    # Names of the regions (also folders names) 
    states = ["CER", "CHA", "CO", "ES", "MO", "SDM", "SP", "ZM"]
    # Varieties in each case
    varieties = [['arabica'], ['arabica'], ['arabica'], ['arabica', 'robusta'], ['arabica'], ['arabica'], ['arabica'], ['arabica', 'robusta']]
            # [['coffee']] 
    # Addresses of the shapefiles with the boundaries for each of the regions
    # Address of the boundary files !!!!SHOULD BE IN EPSG 4326 PROJECTION
    statesBoundFiles = [origin + '/' + s + '/aoi/AOI_' + s + '.shp' for s in states] 
    # Name of subfolder where to save the raw mosaic data (should be in the folders of the regions)
    statesRawFolder = 'raw_data'
    # Name of subfolder where to save the smoothed mosaic data (should be in the folders of the regions)
    statesSmoothFolder = 'smooth_data'
    # Name of subfolder where to save (if produced) and get the reference baseline for each date in the year (should be in the folders of the regions)
    statesRefFolder = 'baseline'
    
    # #DOWNLOAD parameters
    dload = False
    # Product to download
    product = 'MOD13Q1.006'
    # Username for the earthdata website
    user = "olivierpfrancois"
    # Password for the earthdata website
    pwd = "Michele1950"
    # Tiles to download.
    tiles = ['h13v10', 'h13v11', 'h14v10', 'h14v11']  # ['h28v07']
    # Start date for the product download (format YYYY-MM-DD)
    #    If None, will default to date of most recent MODIS file on disk if any, or stop the process
    startDownload = None
    # startDownload = '2017-05-26'
    # End date for the product download (format YYYY-MM-DD)
    #    If None, defaults to today
    endDownload = None
    
    # #MOSAIC parameters
    # Should the downloaded files be mosaiced for each of the regions?
    mosaic = False
    # Starting date for the files to mosaic
    #    If None, will default to the files that have been just downloaded if any.
    startMosaic = '2017-10-01'
    # startMosaic = '2005-01-01'
    # Ending date for the files to mosaic
    #    If None, defaults to today
    endMosaic = None
    # endMosaic = '2005-02-01'
    
    # masking missing values
    checkQuality = False
    inCheck = statesRawFolder
    outCheck = 'masked_missing'
    startCheck = '2011-01-01'
    endCheck = None
    
    # filling missing values
    fillMissing = True
    inMissing = 'masked_missing'
    outMissing = 'filled_missing'
    yearsMissing = [2017]
    daysMissing = [289,305,321,337,353]
    
    # #SMOOTHING parameters
    smooth = False
    regWindow = 7
    avgWindow = 3
    # Starting date for the files to include as input in the smoothing process
    #    If None, defaults to 1 year before the end smoothing date
    startSmooth = None
    # startSmooth = '2012-03-01'
    # Ending date for the files to include as input in the smoothing process
    #    If None, defaults to today
    endSmooth = None
    # Start and end dates for the files to save to the disk after smoothing
    #    If None, defaults to 6 months before end smoothing date
    startSaveS = '2017-06-10'
    # startSaveS = '2017-03-01'
    endSaveS = None  # None to save them up to the end smoothing date
    
    # #Reference rasters parameters
    createBaselines = False
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
    
    # #Ranking individual dates modis images in terms of deciles using the baselines
    ranking = False
    # Starting and ending dates for the images to consider. Included
    #   If None, will default to 60 days before the endRank date
    startRank = '2017-06-15'
    #   If None, will default to today
    endRank = None
    # File to use for masking the output
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
    
    # Minimum density of coffee to consider 
    minCoffee = [0.05, 0.15]
    
    # #Estimate average ndvi value for each region
    avgValue = False
    # Starting and ending dates for the images to consider. Included
    #    If None, defaults to 1 year before today
    startAvg = '2006-01-01'
    #    If None, will default to today
    endAvg = None
    #Grid/Raster to use for the averaging --> FULL PATH!!!!!!!
    # RASTER OPTION
    # avgWeights = [[dst + '/LD/masks/LD_densities_coffee_from_classifications_250m.tif']]
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
    
    weightField = None
    
    #####################################################################################################################
    #####################################################################################################################
    # #ACTIVE CODE
    
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
                print('No start date provided for mosaic. Will mosaic all files on disk')
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
            endCheck= datetime.now()
            endCheck = endCheck.date()
        else:
            endCheck= datetime.strptime(endCheck, '%Y-%m-%d').date()
        
        if startCheck:
            startCheck = datetime.strptime(startCheck, '%Y-%m-%d').date()
        
        for s in states:
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
            
            for d in dataset:
                md.maskQualityVI(ndviRaster=d[0], qualityRaster=d[1], outRaster=d[2], nodataOut=-3000)
    
    if fillMissing:
        
        for s in states:
            inputRasters = [os.path.join(dst, s, inMissing, f) for 
                            f in os.listdir(os.path.join(dst, s, inMissing)) 
                            if f.endswith('.tif')]
            
            inputRasters.sort()
            
            #Get the dates
            datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in inputRasters]
            # Transform into date format
            datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
            # Transform into days from start of the year
            days = [int(d.strftime('%j')) for d in datesAll]
            
            # Get the years for the files on disk
            years = [int(d.strftime('%Y')) for d in datesAll] 
            
            gapfill.gapFill(rasters=inputRasters, seasons=days, years=years,
                            outFolder=os.path.join(dst, s, outMissing),
                            suffix='f', nodata=[-3000], iMax=20,
                            subsetSeasons=daysMissing, subsetYears=yearsMissing, subsetMissing=None,
                            clipRange=(-2000, 10000), parallel=True, nCores=3)
    
    
    if smooth:
        print('Starting the smoothing process')
        md.smoothMODISWrapper(root=dst,
                       regions=states,
                       regionsIn=statesRawFolder,
                       regionsOut=statesSmoothFolder,
                       startSmooth=startSmooth,
                       endSmooth=endSmooth,
                       regWindow=regWindow,
                       avgWindow=avgWindow,
                       startSaveSmooth=startSaveS,
                       endSaveSmooth=endSaveS)
    
    if createBaselines:
        print('Starting creation of baselines')
        md.createBaseline(root=dst, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, regionsOut=statesRefFolder,
                       startRef=startRef, endRef=endRef, mask=maskBaseline, outModelRaster=outModelRaster)
    
    if ranking:
        print('Starting analysis of individual dates compared to baseline')
        md.rankDatesDeciles(root=dst, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, refDecilesIn=statesRefFolder,
                       startRank=startRank, endRank=endRank, mask=maskRank, minDensity=minCoffee)

    if avgValue:
        print('Computing the average ndvi value for each zone')
        md.computeAvgNdvi(root=dst, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, avgWeights=avgWeights,
                       weightField=weightField, startAvg=startAvg, endAvg=endAvg)


if __name__ == '__main__':
    main()
