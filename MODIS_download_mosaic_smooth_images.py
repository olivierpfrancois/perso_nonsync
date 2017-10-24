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

def run_script(iface):
    
    #####################################################################################################################
    #####################################################################################################################
    ##PARAMETERS
    
    #Path to MRT software
    #mrtPath = '/home/olivier/MRT'
    
    ##DIRECTORIES parameters
    #Working directory
    dst = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6' #'E:/gedata_current/jde_coffee/MODIS/collection6'
    #Directory data sources
    origin = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/data' #'E:/gedata_current/jde_coffee/data'
    #Folder to use for temporary files (should be empty)
    tempDir = 'Temp'
    #Destination folder for the download
    rawdataDir = 'raw_data'
    

    ##REGIONS parameters
    #Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    #Names of the regions (also folders names) 
    states = ["CER","CHA","CO","ES","MO","SDM","SP","ZM"]
    #Varieties in each case
    varieties = [['arabica'],['arabica'],['arabica'],['arabica','robusta'],['arabica'],['arabica'],['arabica'],['arabica','robusta']]
    #Addresses of the shapefiles with the boundaries for each of the regions
    #Address of the boundary files !!!!SHOULD BE IN EPSG 4326 PROJECTION
    statesBoundFiles = [origin+'/'+s+'/aoi/AOI_'+s+'.shp' for s in states] 
    #Name of subfolder where to save the raw mosaic data (should be in the folders of the regions)
    statesRawFolder = 'raw_data'
    #Name of subfolder where to save the smoothed mosaic data (should be in the folders of the regions)
    statesSmoothFolder = 'smooth_data'
    #Name of subfolder where to save (if produced) and get the reference baseline for each date in the year (should be in the folders of the regions)
    statesRefFolder = 'baseline'
    
        
    ##DOWNLOAD parameters
    dload = True
    #Product to download
    product = 'MOD13Q1.006'
    #Username for the earthdata website
    user = "olivierpfrancois"
    #Password for the earthdata website
    pwd = "Michele1950"
    #Tiles to download.
    tiles = ['h13v10','h13v11','h14v10','h14v11']
    #Start date for the product download (format YYYY-MM-DD)
    #    If None, will default to date of most recent MODIS file on disk if any, or stop the process
    startDownload = None
    #startDownload = '2017-05-26'
    #End date for the product download (format YYYY-MM-DD)
    #    If None, defaults to today
    endDownload = None
    
    ##MOSAIC parameters
    #Should the downloaded files be mosaiced for each of the regions?
    mosaic = True
    #Starting date for the files to mosaic
    #    If None, will default to the files that have been just downloaded if any.
    startMosaic = None
    #startMosaic = '2005-01-01'
    #Ending date for the files to mosaic
    #    If None, defaults to today
    endMosaic = None
    #endMosaic = '2005-02-01'
    
    
    ##SMOOTHING parameters
    smooth = True
    regWindow = 7
    avgWindow = 3
    #Starting date for the files to include as input in the smoothing process
    #    If None, defaults to 1 year before the end smoothing date
    startSmooth = '2015-04-01'
    #startSmooth = '2016-12-01'
    #Ending date for the files to include as input in the smoothing process
    #    If None, defaults to today
    endSmooth = None
    #Start and end dates for the files to save to the disk after smoothing
    #    If None, defaults to 6 months before end smoothing date
    startSaveS = '2015-09-01'
    #startSaveS = '2017-03-01'
    endSaveS = None #None to save them up to the end smoothing date
    
    
    ##Reference rasters parameters
    createBaselines = False
    #Starting and ending years for the reference period. Included.
    startRef = 2006
    endRef = 2016
    #Mask and output model to use for the baseline rasters. 
    #The mask should overlap perfectly with the modis images
    #The output model can have a different resolution, in which case the baseline will be produced with that resolution
    maskBaseline = [['masks/CER_densities_arabica_from_classifications_250m.tif'],
                    ['masks/CHA_densities_arabica_from_classifications_250m.tif'],
                    ['masks/CO_densities_arabica_from_classifications_250m.tif'],
                    ['masks/ES_densities_arabica_from_classifications_250m.tif','masks/ES_densities_robusta_from_classifications_250m.tif'],
                    ['masks/MO_densities_arabica_from_classifications_250m.tif'],
                    ['masks/SDM_densities_arabica_from_classifications_250m.tif'],
                    ['masks/SP_densities_arabica_from_classifications_250m.tif'],
                    ['masks/ZM_densities_arabica_from_classifications_250m.tif','masks/ZM_densities_robusta_from_classifications_250m.tif']]
    outModelRaster = [['masks/CER_densities_arabica_from_classifications_1km.tif'],
                      ['masks/CHA_densities_arabica_from_classifications_1km.tif'],
                      ['masks/CO_densities_arabica_from_classifications_1km.tif'],
                      ['masks/ES_densities_arabica_from_classifications_1km.tif','masks/ES_densities_robusta_from_classifications_1km.tif'],
                      ['masks/MO_densities_arabica_from_classifications_1km.tif'],
                      ['masks/SDM_densities_arabica_from_classifications_1km.tif'],
                      ['masks/SP_densities_arabica_from_classifications_1km.tif'],
                      ['masks/ZM_densities_arabica_from_classifications_1km.tif','masks/ZM_densities_robusta_from_classifications_1km.tif']]
    
    
    ##Ranking individual dates modis images in terms of deciles using the baselines
    ranking = True
    #Starting and ending dates for the images to consider. Included
    #   If None, will default to 60 days before the endRank date
    startRank = '2017-06-15'
    #   If None, will default to today
    endRank = None
    #File to use for masking the output
    maskRank = [['masks/CER_densities_arabica_from_classifications_1km.tif'],
                ['masks/CHA_densities_arabica_from_classifications_1km.tif'],
                ['masks/CO_densities_arabica_from_classifications_1km.tif'],
                ['masks/ES_densities_arabica_from_classifications_1km.tif','masks/ES_densities_robusta_from_classifications_1km.tif'],
                ['masks/MO_densities_arabica_from_classifications_1km.tif'],
                ['masks/SDM_densities_arabica_from_classifications_1km.tif'],
                ['masks/SP_densities_arabica_from_classifications_1km.tif'],
                ['masks/ZM_densities_arabica_from_classifications_1km.tif','masks/ZM_densities_robusta_from_classifications_1km.tif']]
    #Minimum density of coffee to consider 
    minCoffee = [0.05,0.15]
    
    
    ##Estimate average ndvi value for each region
    avgValue = True
    #Starting and ending dates for the images to consider. Included
    #    If None, defaults to 1 year before today
    startAvg = '2006-01-01'
    #    If None, will default to today
    endAvg = None
    #Grid/Raster to use for the averaging --> FULL PATH!!!!!!!
    '''
    #GRID OPTION
    avgWeights = [origin+'/CER/areas/CER_area_cor.shp',
                  origin+'/CHA/areas/CHA_area_cor.shp',
                  origin+'/CO/areas/CO_area_cor.shp',
                  origin+'/ES/areas/ES_area_cor.shp',
                  origin+'/MO/areas/MO_area_cor.shp',
                  origin+'/SDM/areas/SDM_area_cor.shp',
                  origin+'/SP/areas/SP_area_cor.shp',
                  origin+'/ZM/areas/ZM_area_cor.shp']
    #Name of the variable with the coffee weights if using a shapefile
    #    Put None if using a raster for an area
    weightField = ['arabica_co','arabica_co','arabica_co','robusta_co','arabica_co','arabica_co','arabica_co','arabica_co']
    '''
    #RASTER OPTION
    #avgWeights = [[dst+'/ES/masks/ES_densities_robusta_from_classifications_250m.tif']]
    
    avgWeights = [[dst+'/CER/masks/CER_densities_arabica_from_classifications_250m.tif'],
                  [dst+'/CHA/masks/CHA_densities_arabica_from_classifications_250m.tif'],
                  [dst+'/CO/masks/CO_densities_arabica_from_classifications_250m.tif'],
                  [dst+'/ES/masks/ES_densities_arabica_from_classifications_250m.tif',
                   dst+'/ES/masks/ES_densities_robusta_from_classifications_250m.tif'],
                  [dst+'/MO/masks/MO_densities_arabica_from_classifications_250m.tif'],
                  [dst+'/SDM/masks/SDM_densities_arabica_from_classifications_250m.tif'],
                  [dst+'/SP/masks/SP_densities_arabica_from_classifications_250m.tif'],
                  [dst+'/ZM/masks/ZM_densities_arabica_from_classifications_250m.tif',
                   dst+'/ZM/masks/ZM_densities_robusta_from_classifications_250m.tif']]
    
    weightField = None
    
    #####################################################################################################################
    #####################################################################################################################
    ##ACTIVE CODE
    
    if dload:
        newHdf = downloadMODIS(dstFolder=os.path.join(dst,rawdataDir), 
                               pwd=pwd, user=user, tiles=tiles, product=product, startDownload=startDownload, 
                               endDownload=endDownload)
    else:
        newHdf = []
    
    if mosaic:
        
        thereHdf = [f for f in os.listdir(os.path.join(dst,rawdataDir)) if f.endswith('.hdf')]
        
        if not newHdf:
            if not thereHdf:
                return
            
            if not startMosaic:
                print('No start date to mosaic files already on disk, no mosaic')
        
        if startMosaic:
            startMosaic = datetime.strptime(startMosaic, '%Y-%m-%d').date()
        
        if endMosaic:
            endMosaic = datetime.strptime(endMosaic, '%Y-%m-%d').date()
                
        if newHdf:
            #Dates of these files
            datesNew = [re.search('A([0-9]{7})', f).group(1) for f in newHdf]
            datesNewD = [datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])-1) for d in datesNew]
            if startMosaic:
                startMosaic = min(startMosaic, min(datesNewD))
            else:
                startMosaic = min(datesNewD)
        
        if newHdf or startMosaic:
            print('Starting the mosaic process')
            mosaicMODIS(root=dst, srcFolder=rawdataDir, tmpFolder=tempDir, 
                        regions=states, regionsOut=statesRawFolder, regionsBoundaries=statesBoundFiles, tiles=tiles,
                        startMosaic=startMosaic, endMosaic=endMosaic)
        
    if smooth:
        print('Starting the smoothing process')
        smoothMODIS(root=dst, regions=states, regionsIn=statesRawFolder, regionsOut=statesSmoothFolder, 
                    startSmooth=startSmooth, endSmooth=endSmooth, regWindow=regWindow, avgWindow=avgWindow, 
                    startSaveSmooth=startSaveS, endSaveSmooth=endSaveS)
    
    if createBaselines:
        print('Starting creation of baselines')
        createBaseline(root=dst, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, regionsOut=statesRefFolder, 
                       startRef=startRef, endRef=endRef, mask=maskBaseline, outModelRaster=outModelRaster)
    
    if ranking:
        print('Starting analysis of individual dates compared to baseline')
        rankDatesDeciles(root=dst, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, refDecilesIn=statesRefFolder,
                       startRank=startRank, endRank=endRank, mask=maskRank, minCoffee=minCoffee)
    

    if avgValue:
        print('Computing the average ndvi value for each zone')
        computeAvgNdvi(root=dst, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, avgWeights=avgWeights, 
                       weightField=weightField, startAvg=startAvg, endAvg=endAvg)



#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
##FUNCTIONS


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
        #Defaults to today's date
        endDownload = datetime.now()
        endDownload = endDownload.strftime('%Y-%m-%d')
        
    #Get the names of the hdf files already downloaded and processed
    #thereTif = [f for f in os.listdir(dstFolder) if f.endswith('.tif')]
    thereHdf = [f for f in os.listdir(dstFolder) if f.endswith('.hdf')]
    #Dates of these files
    if thereHdf:
        datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in thereHdf]
        datesHdfD = [datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])) for d in datesHdf]
        #datesTif = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})_', f).group(1) for f in thereTif]
        #datesTif = [datetime.strptime(d, "%Y-%m-%d").date() for d in datesTif]
    else:
        datesHdf = []
        #datesTif = []
    
    #Get latest date of the file on disk
    if datesHdf and not startDownload:
        startDownload = max(datesHdfD).strftime('%Y-%m-%d')
    
    if not startDownload:
        print('There are no files on disk. Choose a start date for the download.')
        return
    
    #Download 
    down = pm.downmodis.downModis(destinationFolder=dstFolder, password=pwd, user=user, 
                url="https://e4ftl01.cr.usgs.gov", tiles=tiles, path='MOLT', product=product, 
                today=endDownload, enddate=startDownload, jpg=False, debug=False, timeout=30, checkgdal=True)
    down._connectHTTP()
    down.downloadsAllDay()
    
    #Get the list of new hdf files on disk
    newHdf = [f for f in os.listdir(dstFolder) if f.endswith('.hdf')]
    newHdf = [f for f in newHdf if f not in thereHdf]
    
    if newHdf:
        print(str(len(newHdf))+' images downloaded')
        
        #Remove the .txt files
        xml = [f for f in dstFolder if f.endswith('.txt') or f.endswith('.log')]
        for f in xml:
            os.remove(os.path.join(dstFolder,f))
    else:
        print('No new images downloaded')
    
    return newHdf


   
def mosaicMODIS(root, srcFolder, tmpFolder, regions, regionsOut, regionsBoundaries, tiles, startMosaic=None, endMosaic=None):
    '''
    Function to mosaic the tiles of modis images and clip them to a series of regions.
    Does not return anything.
    
    root (str): Address of root folder where  srcFolder, tmpFolder, regions folders and regionsOut are located
    srcFolder (str): Address from root where the input downloaded tiles are located
    tmpFolder (str): Address from root of folder where the temporary mosaic before cutting to regions will be stored. 
            !!!!! THIS FOLDER SHOULD BE EMPTY TO START, IT WILL BE EMPTIED AT THE END
    regions (list of str): Names of the regions to process. It should also correspond to the name of the folders
            in root where the regionsOut folders are. 
    regionsOut (str): Name of the folder from regions where the mosaiced and clipped images should be stored. 
            It should be the same for all the regions
    regionsBoundaries (list of str): List the full address of the shapefiles with the regions boundaries
    tiles (list of str): List of  tiles to mosaic (e.g. 'h13v10')
    startMosaic (str): Starting date for the files to mosaic. If None, will process all the files found on the disk.
    endMosaic (str): Ending date for the files to mosaic. If None, defaults to today
    '''
    
    if not endMosaic:
        endMosaic = datetime.now()
        endMosaic = endMosaic.date()
        
    #Get the files on disk
    thereHdf = [f for f in os.listdir(os.path.join(root,srcFolder)) if f.endswith('.hdf')]
    #Dates of these files
    datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in thereHdf]
    datesHdfD = [datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])-1) for d in datesHdf]
    
    #Get the dates to process
    dates = {d for d,D in zip(datesHdf,datesHdfD) if D >= startMosaic and D <= endMosaic}
    dates = sorted(list(dates))
    
    #Get the tiles for the output name
    tilesH = [re.sub('v[0-9]+','',t) for t in tiles]
    tilesH = [re.findall(r'\d+', t) for t in tilesH]
    tilesH = {int(t[0]) for t in tilesH}
    tilesH = list(tilesH)
    tilesH.sort()
    tilesV = [re.sub('h[0-9]+','',t) for t in tiles]
    tilesV = [re.findall(r'\d+', t) for t in tilesV]
    tilesV = {int(t[0]) for t in tilesV}
    tilesV = list(tilesV)
    tilesV.sort()
    
    #Work by date
    for d in dates:
        #Transform into actual date
        dlong = datetime.strptime(d, '%Y%j').date()
        #dlong = datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])-1)
        
        #Progress
        print('Processing date '+dlong.strftime('%Y-%m-%d'))
        
        #Get the file names for that date
        files = [f for f in thereHdf if d in f]
        
        #Check if all the tiles were completed
        complete = True
        for t in tiles:
            if not complete:
                break
            
            complete = any(t in f for f in files)
        
        if not complete:
            print('Not all tiles were downloaded for date '+dlong.strftime('%Y-%m-%d')+". Cannot process.")
            continue
        
        #Create the output name for the mosaic file
        outName = re.sub('\..+','',files[0]) #Get the product name
        outName = outName + '_' + dlong.strftime('%Y-%m-%d') + '_h' + '-'.join(map(str,tilesH)) + 'v' + '-'.join(map(str,tilesV)) + '_250m_16_days_NDVI.tif'
        
        #Mosaic the files and transform to tif
        mos = pm.convertmodis_gdal.createMosaicGDAL(hdfnames=[os.path.join(root,srcFolder,f) for f in files], 
                    subset='1 0 0 0 0 0 0 0 0 0 0 0', outformat='GTiff') #outformat='HDF4Image'
        mos.run(output=os.path.join(root,tmpFolder,'temp.tif')) #'.hdf'
        #mos.write_mosaic_xml(os.path.join(root,tmpFolder,'outName'))
        #mos.write_vrt(os.path.join(root,tmpFolder,outName))
        
        #Change the projection to lat/long
        inRaster = '/'.join([root,tmpFolder,'temp.tif'])
        outRaster = '/'.join([root,tmpFolder,outName])
        cmd = 'gdalwarp -overwrite -t_srs EPSG:4326 -r near -of GTiff -co BIGTIFF=IF_NEEDED %s %s' % (inRaster, outRaster)
        os.system(cmd)
        
        #read the image
        imgMosaic = gdal.Open(outRaster)
        #Get the geotransform and the projection
        geoTrans = imgMosaic.GetGeoTransform()
        projection = imgMosaic.GetProjection()
        #Read band as array
        band = imgMosaic.GetRasterBand(1)
        
        #Loop through the regions to clip and mask the mosaic
        for s, shp in zip(regions, regionsBoundaries):
            
            #Input and output name
            inRaster = root+'/'+tmpFolder+'/'+outName
            outRaster = root+'/'+s+'/'+regionsOut+'/'+outName
            
            # Open the data source and read in the extent
            driver = ogr.GetDriverByName('ESRI Shapefile')
            mask_ds = driver.Open(shp, 0) # 0 means read-only. 1 means writeable.
            mask_layer = mask_ds.GetLayer()
            minX, maxX, minY, maxY = mask_layer.GetExtent()
            
            #Get the pixel position of the shapefile corners
            ulX, ulY = world2Pixel(geoTrans, minX, maxY)
            lrX, lrY = world2Pixel(geoTrans, maxX, minY)
            
            #clip the mosaic image to that extent
            #Import as array xoffset, yoffset, xextent, yextent
            clip = band.ReadAsArray(ulX, ulY, int(lrX - ulX), int(lrY - ulY)).astype(np.float)
            
            #Remove rasters to save if any
            if os.path.isfile(outRaster):
                os.remove(outRaster)
            
            #Create a new geoTransform for the output raster
            geoRegion = list(geoTrans)
            geoRegion[0] = minX
            geoRegion[3] = maxY
            
            #Create an empty raster with the right size
            rows, cols = clip.shape
            driver = gdal.GetDriverByName("GTiff")
            
            new_raster = driver.Create(outRaster, cols, rows, 1, gdal.GDT_Float32, options=['COMPRESS=LZW', 'BIGTIFF=IF_NEEDED'])
            new_raster.SetProjection(projection)
            new_raster.SetGeoTransform(geoRegion)

            new_raster.GetRasterBand(1).SetNoDataValue(-99)
            new_raster.GetRasterBand(1).Fill(-99)
            
            #Rasterize the region into that raster
            gdal.RasterizeLayer(new_raster, [1], mask_layer, burn_values=[1.], options=['ALL_TOUCHED=FALSE'])
            
            #Get the resulting raster band as an array
            new = new_raster.GetRasterBand(1).ReadAsArray().astype(np.float)
            new[new==1.] = clip[new==1.]
            
            #Export the values
            new_raster.GetRasterBand(1).WriteArray(new, 0, 0)
            
            #Close the raster
            new_raster.FlushCache()
            new_raster = None
            
            #Crop the mosaic and mask the raster
            inRaster = root+'/'+tmpFolder+'/'+outName
            outRaster = root+'/'+s+'/'+regionsOut+'/'+outName
            #Call gdalwarp using the os system command
            #cmd = 'gdalwarp -q -cutline %s -crop_to_cutline %s %s' % (shp, inRaster, outRaster)
            #os.system(cmd)
            #Call gdalwarp using the subprocess command
            #out = subprocess.call(['/usr/bin/gdalwarp', '-q', '-overwrite', '-cutline', shp, '-crop_to_cutline', 
            #                       '-of GTIFF', '-co BIGTIFF=IF_NEEDED', inRaster, outRaster])
            #print(out)
            
        imgMosaic = None
        
        #Remove the .xml, .hdf and .vrt intermediary files
        xml = [f for f in os.listdir(os.path.join(root,tmpFolder)) if f.endswith('.xml') or 
                    f.endswith('.hdf') or f.endswith('.vrt') or f.endswith('.tif')]
        for f in xml:
            os.remove(os.path.join(root,tmpFolder,f))



def smoothMODIS(root, regions, regionsIn, regionsOut, startSmooth, endSmooth, regWindow, avgWindow, startSaveSmooth=None, endSaveSmooth=None):
    '''
    Function to do a temporal smoothing of a time series of identical images. Does not return anything, saves smoothed images on disk.
    The input files should contain the date in their name in the format '_([0-9]{4}-[0-9]{2}-[0-9]{2})'
    The processed files will have the same name plus a 'smooth_' prefix
    
    root (str): Address of root folder where  srcFolder, tmpFolder, regions folders and regionsOut are located
    regions (list of str): Names of the regions to process. It should also correspond to the name of the folders
            in root where the regionsOut folders are.
    regionsIn (str): Name of the folder from regions where to find the input files to be smoothed 
            It should be the same for all the regions 
    regionsOut (str): Name of the folder from regions where to save the smoothed images 
            It should be the same for all the regions
    startSmooth (str): Starting date for the files to mosaic. If None, will process all the files found on the disk.
    endSmooth (str): Ending date for the files to mosaic. If None, defaults to today
    regWindow (int): size of the regression window (see Swets et al. for details)
    avgWindow (int): sie of the averaging window (see Swets et al. for details)
    startSaveSmooth (str): Starting date for the files to save. If None, will save all the processed files
    endSaveSmooth (str): Ending date for the files to save. If None, will save all the processed files
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
    
    #Transform into date format
    if not startSaveSmooth:
        startSaveSmooth = endSmooth - timedelta(days=183)
    else:
        startSaveSmooth = datetime.strptime(startSaveSmooth, '%Y-%m-%d').date()
    if not endSaveSmooth:
        endSaveSmooth = endSmooth
    else:
        endSaveSmooth = datetime.strptime(endSaveSmooth, '%Y-%m-%d').date()
    
    #Loop through the regions to do the smoothing for each
    for r in regions:
        print('Processing region '+str(r)+'...')
        
        #Import all the raw modis images on disk
        onDisk = [f for f in os.listdir(os.path.join(root,r,regionsIn)) if f.endswith('.tif')]
        
        #Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        #Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        
        #Keep only the files and dates within the dates to process
        onDisk = [f for f,d in zip(onDisk,datesAll) if d >= startSmooth and d <= endSmooth]
        datesAll = [d for d in datesAll if d >= startSmooth and d <= endSmooth]
        
        #Sort the two list by date
        datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
        #Create a mask to know which of the dates to save
        toSave = [d >= startSaveSmooth and d <= endSaveSmooth for d in datesAll]
        
        #Import all the images to process in the right order
        toProcess = [gdal.Open(os.path.join(root,r,regionsIn,f)) for f in onDisk]
        
        #Get the no data value if any
        nodata = toProcess[1].GetRasterBand(1).GetNoDataValue()
        
        #Remove rasters to save if any
        for s, f in zip(toSave, onDisk):
            if s and os.path.isfile(root+'/'+r+'/'+regionsOut+'/smooth_'+f):
                os.remove(root+'/'+r+'/'+regionsOut+'/smooth_'+f)
        
        #Create empty copies of these files to use for the smoothed data only for the files to save
        processed = [new_raster_from_base(p, root+'/'+r+'/'+regionsOut+'/smooth_'+f, 'GTiff', np.nan, gdal.GDT_Float32) 
                        if s else False for s,p,f in zip(toSave, toProcess, onDisk)]
        
        #Get the size of the rasters to identify the limits of the blocks to loop through
        band = toProcess[1].GetRasterBand(1)
        #Get the size of the raster
        xsize = band.XSize
        ysize = band.YSize
        #Set the block size
        BlockXSize = 256
        BlockYSize = 256
        #Get the number of blocks in x and y directions
        xBlocks = int(round(xsize/BlockXSize)) + 1
        yBlocks = int(round(ysize/BlockYSize)) + 1
        
        totBlocks = xBlocks*yBlocks
        progress = 1
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                print('   '+str(r)+': Processing block '+str(progress)+' of '+str(totBlocks))
                progress += 1
                
                blocks = [readRasterBlock(p, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize) for p in toProcess]
                
                #Bring the blocks together into one single array
                blocks = np.dstack(blocks)
                #Recast the type
                blocks = blocks.astype(np.float32)
                
                #Do the smoothing for each pixel
                blocks = smoothingSwets(blocks, regWindow, avgWindow, nodata)
                #blocks = smoothingSavitzky(blocks, nodata)
                blocks = np.divide(blocks,10000.)
                
                #Change the values in the output raster
                blocks = np.dsplit(blocks, len(toProcess))
                for s, p, b in zip(toSave, processed, blocks):
                    if s:
                        p.GetRasterBand(1).WriteArray(b[:,:,0], xStep*BlockXSize, yStep*BlockYSize)
        
        #Close the rasters
        for s, p in zip(toSave, processed):
            if s:
                p.FlushCache()
                p = None
        
        for p in zip(toSave, toProcess):            
            p=None



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
            pixel = np.copy(block[X,Y,:])
            
            #Check if entire pixel vector is nan, just return it
            if nodata and np.all(pixel==nodata):
                pixel.fill(np.nan)
                block[X,Y,:] = pixel
                continue
            elif np.all(np.isnan(pixel)):
                continue
            
            #Get the shape of the original data
            originShape = pixel.shape
            
            #Reshape the data
            pixel = np.reshape(pixel, (len(pixel),1))
            pixel = pixel.astype(np.float32)
            
            #Replace no data values by nan if needed
            if nodata:
                pixel[pixel==nodata] = np.nan
            
            #Interpolate missing values if any
            #THIS SHOULD BE DONE BEFOREHAND FOR THE ENTIRE RASTER USING TEMPORAL AND SPATIAL INFERENCE
            #Get a mask of the nan values
            if np.isnan(pixel).any():
                pixel = np.reshape(pixel, len(pixel))
                nans = np.isnan(pixel)
                #Get the actual indices for these values and the non nan values
                notNans = ~nans
                notNans = notNans.nonzero()[0]
                nans = nans.nonzero()[0]
                #Interpolate
                pixel[nans]= np.interp(nans, notNans, pixel[notNans])
                #Reshape
                pixel = np.reshape(pixel, (len(pixel),1))
            
            #Create a vector of weights for the regression
            weights = np.copy(pixel)
            weights[0,0] = 0.5 #Assume first value is a middle point
            for i in range(1, len(pixel)-1, 1):
                #Local high values get a weight of 1.5, local middle values get a weight of 0.5 and local low values only 0.005
                if pixel[i-1,0] < pixel[i,0] and pixel[i,0] > pixel[i+1,0]:
                    weights[i,0] = 1.5
                elif ((pixel[i-1,0] <= pixel[i,0] and pixel[i,0] <= pixel[i+1,0]) or (pixel[i-1,0] >= pixel[i,0] and pixel[i,0] >= pixel[i+1,0])):
                    weights[i,0] = 0.5
                elif pixel[i-1,0] > pixel[i,0] and pixel[i,0] < pixel[i+1,0]:
                    weights[i,0] = 0.005
            #For the last point
            if pixel[len(pixel)-1,0] >= pixel[len(pixel)-2,0] and pixel[len(pixel)-1,0] >= pixel[len(pixel)-3,0]:
                #If the last data point is greater than the previous 2, then assume it is a high point
                weights[len(pixel)-1,0] = 1.5
            elif pixel[len(pixel)-1,0] >= pixel[len(pixel)-2,0] and pixel[len(pixel)-1,0] < pixel[len(pixel)-3,0]:
                #If the last data point is greater than the previous only, then assume it is a middle point
                weights[len(pixel)-1,0] = 0.5
            elif (pixel[len(pixel)-1,0] < pixel[len(pixel)-2,0] and pixel[len(pixel)-1,0] >= pixel[len(pixel)-3,0] and 
                  weights[len(pixel)-3,0] != 0.005):
                #If less than the previous but more than the one before that is not a low point
                weights[len(pixel)-1,0] = 0.5
            else:
                #If the last point is less than the last 2, or greater than only one of the two, assume it is a low point
                weights[len(pixel)-1,0] = 0.002
            
            #Create a matrix with the data for this pixel and for the weights
            #For the data:
            #Each column will be the same
            dataRaw = np.repeat(pixel, len(pixel)+regWindow, axis=1)
            #Set to nan all data for each column besides the data used for the regressions
            ltri = np.tril_indices(n=len(pixel), k=-1, m=len(pixel)+regWindow) #Lower triangle indices below the diagonal
            dataRaw[ltri] = np.nan
            utri = np.triu_indices(n=len(pixel), k=0, m=len(pixel))
            dataRaw[:,(regWindow):][utri] = np.nan
            #Remove the first two and last 3 columns, since they don't have enough points for a regression
            dataRaw = dataRaw[:,2:(len(pixel)+regWindow-3)]
            
            #For the weights:
            weights = np.repeat(weights, len(pixel)+regWindow-5, axis=1)
            weights[np.isnan(dataRaw)] = np.nan
            
            #Create an empty array for the results of the regressions
            dataReg = np.zeros(dataRaw.shape)
            
            #Estimate the regression for each column of dataRaw
            for i in range(len(pixel)+regWindow-5):
                #Prepare regression data
                y = dataRaw[:,i][~np.isnan(dataRaw[:,i])] #dependent
                w = weights[:,i][~np.isnan(dataRaw[:,i])] #weights
                x = range(1,len(y)+1) #independent
                x = np.asarray(x)
                x = x.astype(np.float32)
                x = np.reshape(x, y.shape)
                #Estimate potential outliers
                out = np.divide(np.absolute(y-np.mean(y)),np.std(y)) #Will be an outlier if greater than 3 (only extreme outliers should be picked up)
                out = np.reshape(out, y.shape)
                #Remove outliers before regression
                yout = y[out<3]
                wout = w[out<3]
                xout = x[out<3]
                #Compute parameters of regression
                numerator = np.sum(wout)*np.sum(np.multiply(np.multiply(wout,xout),yout))-np.sum(np.multiply(wout,yout))*np.sum(np.multiply(wout,xout))
                denominator = np.sum(wout)*np.sum(np.multiply(wout,np.square(xout)))-np.square(np.sum(np.multiply(wout,xout)))
                b = np.divide(numerator, denominator)
                numerator = (np.sum(np.multiply(wout,yout))-b*np.sum(np.multiply(wout,xout)))
                denominator = np.sum(wout)
                a = np.divide(numerator, denominator)
                #Compute the predicted values from the regression
                dataReg[:,i][~np.isnan(dataRaw[:,i])] = a + b*x
            
            dataReg[np.isnan(dataRaw)] = np.nan
            
            #Combination of the results from the regression for each point
            #Now we combine for each point the results from the moving regression window
            #We take into account the results from avg.window regression windows, centered around the point, unless we are at the edges
            smoothed = [] #Will hold the smoothed results
            t = int(np.floor(avgWindow/2.)) #number of predicted values from regressions to take into account around the center
            for i in range(len(pixel)):
                x = dataReg[i,:][~np.isnan(dataReg[i,:])]
                if i < np.floor(regWindow/2.) and len(x) < regWindow:
                    center = int(np.floor(regWindow/2.))-(regWindow-len(x))
                    res = np.mean(x[max(0,center-t):max(1,center+t+1)])
                elif i == len(pixel)-1:
                    res = np.mean(x[(int(np.floor(regWindow/2.))-t-1):(int(np.floor(regWindow/2.))+t+1)])
                else:
                    res = np.mean(x[(int(np.floor(regWindow/2.))-t):(int(np.floor(regWindow/2.))+t+1)])    
                smoothed.append(res)
            
            smoothed = np.asarray(smoothed)
            smoothed = np.reshape(smoothed, originShape)
            
            block[X,Y,:] = smoothed
    
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
            pixel = np.copy(block[X,Y,:])
            
            #Check if entire pixel vector is nan, just return it
            if nodata and np.all(pixel==nodata):
                pixel.fill(np.nan)
                block[X,Y,:] = pixel
                continue
            elif np.all(np.isnan(pixel)):
                continue
            
            #Get the shape of the original data
            originShape = pixel.shape
            
            #Reshape the data
            pixel = np.reshape(pixel, len(pixel))
            pixel = pixel.astype(np.float32)
            
            #Replace no data values by nan if needed
            if nodata:
                pixel[pixel==nodata] = np.nan
            
            #Get a copy of the pixel with only the non nan values
            toSmooth = pixel[~np.isnan(pixel)]
            
            #The moving window of 7 doesn't work if there aren't enough data
            if len(toSmooth) < 7:
                pixel.fill(np.nan)
                block[X,Y,:] = pixel
                continue
            
            for _ in range(3):
                #Take a first smoothing of the pixel
                smoothed = savitzkyAvg(toSmooth)
                #Replace the values that are lower in the original pixel by the values in the original smoothed pixel
                np.copyto(toSmooth, smoothed, where=toSmooth<smoothed)
            
            #Replace the smoothed values in the original pixel
            pixel[~np.isnan(pixel)] = toSmooth
            
            #Interpolate missing values if any
            #Get a mask of the nan values
            if np.isnan(pixel).any():
                nans = np.isnan(pixel)
                #Get the actual indices for these values and the non nan values
                notNans = ~nans
                notNans = notNans.nonzero()[0]
                nans = nans.nonzero()[0]
                #Interpolate
                pixel[nans]= np.interp(nans, notNans, pixel[notNans])
                
            #Reshape
            pixel = np.reshape(pixel, originShape)
            
            #Replace the smoothed pixel in the block
            block[X,Y,:] = pixel
    
    return block

def savitzkyAvg(pixel):
    #This function averages a pixel using a window 7 for the Savitzky Golay and returns a pixel of the same shape
    
    out = np.copy(pixel)
    
    for i in range(len(pixel)):
        if i < 3:
            #Take the moving window
            avg = pixel[0:7]
            #Do the weighted sum
            if i == 0:
                out[i] = np.sum(np.multiply(avg, np.array([0.1190,-0.0714,-0.1429,-0.0952,0.0714,0.3571,0.7619])))
            elif i == 1:
                out[i] = np.sum(np.multiply(avg, np.array([-0.0714,-0.0000,0.0714,0.1429,0.2143,0.2857,0.3571])))
            elif i == 2:
                out[i] = np.sum(np.multiply(avg, np.array([-0.1429,0.0714,0.2143,0.2857,0.2857,0.2143,0.0714])))
        elif i > (len(pixel)-4):
            #Take the moving window
            avg = pixel[(len(pixel)-7):]
            #Do the weighted sum
            if i == len(pixel)-3:
                out[i] = np.sum(np.multiply(avg, np.array([0.0714,0.2143,0.2857,0.2857,0.2143,0.0714,-0.1429])))
            elif i == len(pixel)-2:
                out[i] = np.sum(np.multiply(avg, np.array([0.3571,0.2857,0.2143,0.1429,0.0714,-0.0000,-0.0714])))
            elif i == len(pixel)-1:
                out[i] = np.sum(np.multiply(avg, np.array([0.7619,0.3571,0.0714,-0.0952,-0.1429,-0.0714,0.1190])))
        else:
            #Take the moving window
            avg = pixel[i-3:i+4]
            #Do the weighted sum
            out[i] = np.sum(np.multiply(avg, np.array([-0.0952,0.1429,0.2857,0.3333,0.2857,0.1429,-0.0952])))                
        
    return(out)
            

def createBaseline(root, regions, varieties, regionsIn, regionsOut, startRef, endRef, mask=None, outModelRaster=None):
    #Loop through the regions
    for r in range(len(regions)):
        #Get the names of all the smoothed rasters on file
        onDisk = [f for f in os.listdir(os.path.join(root,regions[r],regionsIn)) if f.endswith('.tif')]
        
        #Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        #Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        
        #Transform into days from start of the year and keep only the unique values
        days = {d.strftime('%j') for d in datesAll}
        #Transform back into a list and order by date
        days = [int(d) for d in days]
        days.sort()
        
        for v in range(len(varieties[r])):
            #Loop through the dates to create a baseline raster for each
            for d in days:
                #Get the names of all the rasters for this date
                dates = [datetime(y, 1, 1).date()+timedelta(days=d-1) for y in range(startRef,endRef+1,1)]
                #dates = [datetime.strptime('0101'+str(y), '%d%m%Y').date()+timedelta(days=d-1) for y in range(startRef,endRef+1,1)]
                files = [f for f,date in zip(onDisk,datesAll) if date in dates]
                #Import all the images to use for estimating the deciles
                if outModelRaster and outModelRaster[r]:
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
                    if mask and mask[r]:
                        #I need to mask the rasters first to only have the coffee pixels when I change the resolution
                        #Import the mask raster as an array
                        p = gdal.Open(os.path.join(root,regions[r],mask[r][v]))
                        nanMask = p.GetRasterBand(1).ReadAsArray()
                        #Transform all the non zero values to nan
                        nanMask = np.logical_or(nanMask == 0, np.isnan(nanMask))
                    
                    #Mask the rasters
                    toProcess = []
                    for f in files:
                        p = gdal.Open(os.path.join(root,regions[r],regionsIn,f))
                        pBand = p.GetRasterBand(1).ReadAsArray()
                        if mask and mask[r]:
                            pBand[nanMask] = np.nan
                            
                        toProcess.append(new_raster_from_base(p, '', 'MEM', np.nan, gdal.GDT_Float32, bands=1))
                        toProcess[-1].GetRasterBand(1).WriteArray(pBand)
                    p = None
                    pBand = None
                    
                    #Change the resolution of the rasters
                    outModel = gdal.Open(os.path.join(root,regions[r],outModelRaster[r][v]))
                    toProcess = [warp_raster(p, outModel, resampleOption='average', outputURI=None, outFormat='MEM') 
                                 for p in toProcess]
                else:
                    #Import the images
                    toProcess = [gdal.Open(os.path.join(root,regions[r],regionsIn,f)) for f in files]
                    
                #Get the no data value if any
                nodata = toProcess[0].GetRasterBand(1).GetNoDataValue()
                
                #Name of the output baseline for that date
                outname = 'ndvi_deciles_0to100pct_ref_period-'+str(startRef)+'-'+str(endRef)+'_day'+str(d)+'_date-'+dates[0].strftime('%b-%d')+'_'+varieties[r][v]+'.tif'
                
                #Remove existing raster if any
                if os.path.isfile(root+'/'+regions[r]+'/'+regionsOut+'/'+outname):
                    os.remove(root+'/'+regions[r]+'/'+regionsOut+'/'+outname)
                
                #Create an empty copy with 10 layers to use for storing the deciles
                processed = new_raster_from_base(toProcess[0], root+'/'+regions[r]+'/'+regionsOut+'/'+outname, 'GTiff', np.nan, gdal.GDT_Float32, bands=10)
                
                #Get the size of the rasters to identify the limits of the blocks to loop through
                band = toProcess[1].GetRasterBand(1)
                #Get the size of the raster
                xsize = band.XSize
                ysize = band.YSize
                #Set the block size
                BlockXSize = 256
                BlockYSize = 256
                #Get the number of blocks in x and y directions
                xBlocks = int(round(xsize/BlockXSize)) + 1
                yBlocks = int(round(ysize/BlockYSize)) + 1
                
                for xStep in range(xBlocks):
                    for yStep in range(yBlocks):
                    
                        block = [readRasterBlock(p, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize) for p in toProcess]
                        
                        #Bring the blocks together into one single array
                        block = np.dstack(block)
                        #Recast the type
                        block = block.astype(np.float32)
                        
                        #Estimate the deciles for each pixel
                        deciles = estimateDeciles(block, nodata)
                        
                        #Change the values in the output raster
                        deciles = np.dsplit(deciles, 10)
                        for i in range(10):
                            processed.GetRasterBand(i+1).WriteArray(deciles[i][:,:,0], xStep*BlockXSize, yStep*BlockYSize)
                
                #Close the rasters
                for p in toProcess:
                    p.FlushCache()
                    p = None
                processed.FlushCache()
                processed = None
    
def estimateDeciles(block, nodata):
    extent = block.shape
    deciles = np.empty((extent[0],extent[1],10))

    for X in range(extent[0]):
        for Y in range(extent[1]):
            pixel = np.copy(block[X,Y,:])
            
            #Replace the no data value by nan if nodata was provided
            if nodata:
                pixel[pixel==nodata] = np.nan
            
            #Get the number of nan values in the pixel
            nbNodata = np.sum(np.isnan(pixel))
            
            #Return nan if there are not more than 6 valid values
            if len(pixel) - nbNodata < 6:
                deciles[X,Y,:].fill(np.nan)
                continue
            
            #Remove the nan values and reshape the data
            pixel = pixel[~np.isnan(pixel)]
            pixel = np.reshape(pixel, (len(pixel),1))
            pixel = pixel.astype(np.float32)
            
            #Find the deciles
            deciles[X,Y,:] = np.reshape(np.percentile(pixel, np.arange(0, 100, 10)),deciles[X,Y,:].shape)
    
    return deciles



def rankDatesDeciles(root, regions, varieties, regionsIn, refDecilesIn, startRank, endRank, mask, minCoffee):
    
    #Transform into date format
    if not endRank:
        endRank = datetime.now()
        endRank = endRank.date()
    else:
        endRank = datetime.strptime(endRank, '%Y-%m-%d').date()
    if not startRank:
        startRank = endRank - timedelta(days=60)
    else:
        startRank = datetime.strptime(startRank, '%Y-%m-%d').date()
    
    #Loop through the regions to do the ranking for each
    for r in range(len(regions)):
        print('Ranking region '+str(regions[r])+'...')
        
        #Import all the smooth modis images on disk
        onDisk = [f for f in os.listdir(os.path.join(root,regions[r],regionsIn)) if f.endswith('.tif')]
        
        #Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        #Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
        #Keep only the files and dates within the dates to process
        onDisk = [f for f,d in zip(onDisk,datesAll) if d >= startRank and d <= endRank]
        datesAll = [d for d in datesAll if d >= startRank and d <= endRank]
        #Sort the two list by date
        datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
        #Transform into days from start of the year and keep only the unique values
        days = [int(d.strftime('%j')) for d in datesAll]
        
        #Get the names of the decile baselines on disk
        baseFiles = [f for f in os.listdir(os.path.join(root,regions[r],refDecilesIn)) if f.endswith('.tif')]
        
        for v in range(len(varieties[r])):
            #Import the mask
            nanMask = gdal.Open(os.path.join(root,regions[r],mask[r][v]))
            
            #Loop through the days to rank each date
            for day, date, fileS in zip(days,datesAll,onDisk):
                #Import the smoothed image
                smoothImg = gdal.Open(os.path.join(root,regions[r],regionsIn,fileS))
                
                nodata = smoothImg.GetRasterBand(1).GetNoDataValue()
                
                #Import the decile baseline
                baseFile = [f for f in baseFiles if 'day'+str(day)+'_' and varieties[r][v] in f]
                baseImg = gdal.Open(os.path.join(root,regions[r],refDecilesIn,baseFile[0]))
                
                #Adapt the resolution of the smooth images to the baseline if needed
                geoSmooth = smoothImg.GetGeoTransform()
                geoBase = baseImg.GetGeoTransform()
                reproject = [1 for a,b in zip(geoSmooth,geoBase) if not a==b]
                if reproject:
                    smoothImg = warp_raster(smoothImg, baseImg, resampleOption='average', outputURI=None, outFormat='MEM')
                
                #Get the size of the rasters to identify the limits of the blocks to loop through
                band = baseImg.GetRasterBand(1)
                #Get the size of the raster
                xsize = band.XSize
                ysize = band.YSize
                #Set the block size
                BlockXSize = 256
                BlockYSize = 256
                #Get the number of blocks in x and y directions
                xBlocks = int(round(xsize/BlockXSize)) + 1
                yBlocks = int(round(ysize/BlockYSize)) + 1
                band=None
                
                for m in minCoffee:
                    #Create the output name
                    outname = 'ndvi_'+date.strftime('%Y-%m-%d')+'_CompareToDecile_0BelowMin_110AboveMax_maskedbelow'+str(int(m*100))+'%'+varieties[r][v]+'.tif'
                    #Remove existing raster if any
                    if os.path.isfile(root+'/'+regions[r]+'/'+outname):
                        os.remove(root+'/'+regions[r]+'/'+outname)
                    #Create an empty copy for storing the deciles comparisons
                    processed = new_raster_from_base(baseImg, root+'/'+regions[r]+'/'+outname, 'GTiff', -32768, gdal.GDT_Int16, bands=1)
                    
                    for xStep in range(xBlocks):
                        for yStep in range(yBlocks):
                            
                            #Read the block from the images
                            blockSmooth = readRasterBlock(smoothImg, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize)
                            
                            blockBase = [readRasterBlock(baseImg, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize, band=b+1) 
                                         for b in range(baseImg.RasterCount)]
                            
                            #Read the block from the mask 
                            blockMask = readRasterBlock(nanMask, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize)
                            blockMask = blockMask < m
                            
                            #Bring the blocks together into one single array
                            blockBase = np.dstack(blockBase)
                            
                            #Apply the mask
                            blockSmooth[blockMask] = np.nan
                            blockBase[blockMask] = np.nan
                            
                            #Recast the type to be sure
                            blockSmooth = blockSmooth.astype(np.float32)
                            blockBase = blockBase.astype(np.float32)
                        
                            #Estimate the placement for each pixel
                            ranks = estimateRank(block=blockSmooth, ref=blockBase, nodata=nodata)
                            
                            #Change the values in the output raster
                            processed.GetRasterBand(1).WriteArray(ranks, xStep*BlockXSize, yStep*BlockYSize)
                
                #Close the rasters
                smoothImg.FlushCache()
                smoothImg = None
                baseImg.FlushCache()
                baseImg = None
                processed.FlushCache()
                processed = None

def estimateRank(block, ref, nodata):
    extent = block.shape
    ranking = np.empty((extent[0],extent[1]), dtype=int)
    
    for X in range(extent[0]):
        for Y in range(extent[1]):
            testPixel = np.copy(block[X,Y])
            refPixel = np.copy(ref[X,Y,:])
            refPixel = np.reshape(refPixel, (len(refPixel),1))
            
            #Replace the no data value by nan if nodata was provided
            if nodata:
                refPixel[np.logical_or(refPixel==nodata,np.isinf(refPixel))] = np.nan
                if testPixel==nodata or np.isinf(testPixel):
                    testPixel = np.nan
            else:
                refPixel[np.isinf(refPixel)] = np.nan
                if np.isinf(testPixel):
                    testPixel = np.nan
            
            #return nan if the test pixel is nodata
            if np.isnan(testPixel) or np.all(np.isnan(refPixel)):
                ranking[X,Y] = -32768
                continue
            
            #Get the position of the pixel in the reference deciles
            refPixel = np.reshape(refPixel, len(refPixel))
            
            rank = np.searchsorted(refPixel, testPixel)*10
            
            #If the value is above the deciles, returns length of the vector
            if len(refPixel)==10 and rank==100:
                rank = 110
            #If the value is at min or below, returns 0
            if rank==0:
                if testPixel == refPixel[0]:
                    rank = 10
                elif testPixel/refPixel[0]>=0.75:
                    rank = -1
                elif testPixel/refPixel[0]>=0.5:
                    rank=-2
                else:
                    rank=-3
            
            ranking[X,Y] = rank
    
    return np.reshape(ranking, extent)



def computeAvgNdvi(root, regions, varieties, regionsIn, avgWeights, weightField=None, startAvg=None, endAvg=None, alltouch=False):
    
    #Transform into date format
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
    
    #Create an empty dictionary to get the values for each of the regions
    averages = {}
    colnames = []
    
    for r in range(len(regions)):
        for v in range(len(varieties[r])):
            if not avgWeights[r][v]:
                continue
            
            colnames.append(regions[r]+'_'+varieties[r][v])
            
            #Get the images to consider
            print('Averaging region '+str(regions[r])+'...')
            
            #Import all the smooth modis images on disk
            onDisk = [f for f in os.listdir(os.path.join(root,regions[r],regionsIn)) if f.endswith('.tif')]
            
            #Dates of these files
            datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
            #Transform into date format
            datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
            #Keep only the files and dates within the dates to process
            onDisk = [f for f,d in zip(onDisk,datesAll) if d >= startAvg and d <= endAvg]
            datesAll = [d for d in datesAll if d >= startAvg and d <= endAvg]
            
            if not onDisk:
                print('no modis images to process in '+regions[r])
                continue
            
            #Get a base image as a reference for format
            baseImg = gdal.Open(os.path.join(root,regions[r],regionsIn,onDisk[0]))
            
            #Rasterize the gridded weights if needed or simply import it
            if avgWeights[r][v].endswith('.shp'):
                if not weightField[r][v]:
                    print('The name of the field with the densities needs to be specified to average')
                    break
                
                #Import the vector layer
                #Open the shapefile
                driver = ogr.GetDriverByName('ESRI Shapefile')
                dataSource = driver.Open(avgWeights[r][v], 0) # 0 means read-only. 1 means writeable.
            
                # Create layer
                inVecLayer = dataSource.GetLayer(0)
                                        
                # Prepare an empty raster to rasterize the shapefile
                weightsRaster = new_raster_from_base(baseImg, 'temp', 'MEM', -1, gdal.GDT_Float32) 
                
                #Transform alltouch
                if alltouch:
                    alltouch = 'TRUE'
                else:
                    alltouch = 'FALSE'
            
                # Rasterize the vector layer:
                gdal.RasterizeLayer(weightsRaster, [1], inVecLayer, options=['ALL_TOUCHED='+alltouch, 'ATTRIBUTE='+weightField[r][v]])
                
                inVecLayer = None
                
            elif avgWeights[r][v].endswith(('.tif','.TIF')):
                weightsRaster = gdal.Open(avgWeights[r][v])
                
                #Change the resolution of the raster to match the images if needed
                
                geoSmooth = weightsRaster.GetGeoTransform()
                geoBase = baseImg.GetGeoTransform()
                reproject = [1 for a,b in zip(geoSmooth,geoBase) if not a==b]
                if reproject:
                    weightsRaster = warp_raster(weightsRaster, baseImg, resampleOption='nearest', outputURI=None, outFormat='MEM')
            
            baseImg = None
            nodataWeights = weightsRaster.GetRasterBand(1).GetNoDataValue()
            
            #Loop through the images to compute the average for each
            for img, date in zip(onDisk,datesAll):
                
                #Import the image
                baseImg = gdal.Open(os.path.join(root,regions[r],regionsIn,img))
                
                #Loop through the blocks to compute the average for each
                
                #Get the size of the rasters to identify the limits of the blocks to loop through
                band = baseImg.GetRasterBand(1)
                #Get the no data value
                nodataBase = band.GetNoDataValue()
                #Get the size of the raster
                xsize = band.XSize
                ysize = band.YSize
                #Set the block size
                BlockXSize = 256
                BlockYSize = 256
                #Get the number of blocks in x and y directions
                xBlocks = int(round(xsize/BlockXSize)) + 1
                yBlocks = int(round(ysize/BlockYSize)) + 1
                band=None
                
                sumNdvi = []
                sumWeights = []
                
                for xStep in range(xBlocks):
                    for yStep in range(yBlocks):
                        
                        #Read the block from the image
                        blockBase = readRasterBlock(baseImg, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize)
                        
                        #Read the block from the weights
                        blockWeight = readRasterBlock(weightsRaster, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize)
                        
                        #Recast the type to be sure
                        blockBase = blockBase.astype(np.float32)
                        blockWeight = blockWeight.astype(np.float32)
                        
                        #Replace the no data values by 0
                        blockBase[np.logical_or(np.logical_or(blockBase == nodataBase,np.isnan(blockBase)),np.logical_or(blockBase>1,blockBase<-1))] = 0.
                        blockWeight[np.logical_or(np.logical_or(blockWeight == nodataWeights,np.isnan(blockWeight)),np.logical_or(blockWeight>1,blockWeight<0))] = 0.
                    
                        #Estimate the weighted sum for each pixel
                        sumNdvi.append(np.sum(np.multiply(blockBase,blockWeight)))
                        sumWeights.append(np.sum(blockWeight))
                        
                #Combine for the entire image
                sumNdvi = np.sum(sumNdvi)/np.sum(sumWeights)
                #Transform into per Hectare
                sumNdvi = sumNdvi/(250)*10000 
                
                #Add to the output dictionary along with the date and region as a dictionary
                #averages.append({'region':regions[r],'date':date,'ndviPerHa':sumNdvi})
                if date.strftime('%Y-%m-%d') in averages:
                    averages[date.strftime('%Y-%m-%d')][regions[r]+'_'+varieties[r][v]] = sumNdvi
                else:
                    averages[date.strftime('%Y-%m-%d')] = {}
                    averages[date.strftime('%Y-%m-%d')]['date'] = date.strftime('%Y-%m-%d')
                    averages[date.strftime('%Y-%m-%d')][regions[r]+'_'+varieties[r][v]] = sumNdvi
            
    #Export the dictionary
    outNm = 'Weighted_avg_ndvi_'+startAvg.strftime('%Y-%m-%d')+'_'+endAvg.strftime('%Y-%m-%d')+'.txt'
    #Sort the dates
    datesAll.sort()
    #order the output by date in a list. Each element is an element of the original dictionary and will be exported
    out = []
    for date in datesAll:
        out.append(averages[date.strftime('%Y-%m-%d')])
    with open(os.path.join(root,outNm), "w") as f:
        dict_writer = DictWriter(f, ['date']+colnames, extrasaction='ignore', delimiter="\t", restval="0")
        dict_writer.writeheader()
        for p in out:
            dict_writer.writerow(p)


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
    
    #Define resampling options
    resampleOptions = {'nearest': gdalconst.GRA_NearestNeighbour, 'bilinear':gdalconst.GRA_Bilinear, 
                   'cubic':gdalconst.GRA_Cubic, 'cubic spline':gdalconst.GRA_CubicSpline, 
                   'lanczos':gdalconst.GRA_Lanczos, 'average':gdalconst.GRA_Average, 
                   'mode':gdalconst.GRA_Mode} 
    
    if not resampleOption in resampleOptions.keys():
        return False
    
    #Raster to host the warped output
    if outFormat == 'MEM':
        rOut = new_raster_from_base(dst, 'temp', 'MEM', 0, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    else:
        rOut = new_raster_from_base(dst, outputURI, outFormat, 0, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    
    #Warp: the parameters are source raster, destination raster, source projection, destination projection, resampling option 
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
    
    #Get the wanted band
    band = src.GetRasterBand(band)
    
    #Get the size of the raster
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
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line)