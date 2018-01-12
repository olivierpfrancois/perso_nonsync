# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Script to download MODIS images from the server, 
    and mosaic the images and mask them for each region
    to process.
 """

import os, re
import MODIS_gedata_toolbox as md


def main():
    #####################################################################################################################
    #####################################################################################################################
    # #PARAMETERS
    
    # #DIRECTORIES parameters
    # Working directory
    dst = '/home/olivierp/jde_coffee/MODIS/collection6/Brazil'
    # /media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6' #'E:/gedata_current/jde_coffee/MODIS/collection6'
    # Directory data sources
    origin = '/home/olivierp/jde_coffee/data/Brazil'
    # '/media/olivier/olivier_ext1/gedata_current/jde_coffee/data' #'E:/gedata_current/jde_coffee/data'
    # Folder inside dst to use for temporary files (should be empty)
    tempDir = '/home/olivierp/jde_coffee/Temp'
    # Destination folder for the download
    rawdataDir = 'raw_data'
    
    # #REGIONS parameters
    # Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    # Names of the regions (also folders names) 
    states = ["CER", "CHA", "CO", "ES", "MO", "SDM", "SP", "ZM"]  # ["LD"]
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
    startMosaic = '2005-01-01'
    # startMosaic = '2005-01-01'
    # Ending date for the files to mosaic
    #    If None, defaults to today
    endMosaic = None
    # endMosaic = '2005-02-01'
    
    # masking missing values
    checkQuality = True
    inFolder = statesRawFolder
    outFolder = 'masked_missing'
    
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
                                  startMosaic=startMosaic,
                                  endMosaic=endMosaic)
    
    if checkQuality:
        for s in states:
            # Get all the images on disk
            allNDVI = [os.path.join(dst, s, inFolder, f) for 
                       f in os.listdir(os.path.join(dst, s, inFolder)) 
                       if f.endswith('NDVI.tif')]
            allNDVI.sort()
            allQuality = [os.path.join(dst, s, inFolder, f) for 
                          f in os.listdir(os.path.join(dst, s, inFolder)) 
                          if f.endswith('Quality.tif')]
            allQuality.sort()
            
            allOut = [re.sub(inFolder, outFolder, f) for f in allNDVI]
            
            # Define the dataset
            dataset = zip(allNDVI, allQuality, allOut)
            
            for d in dataset:
                md.maskQualityVI(d[0], d[1], d[2])
    
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
