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

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
import pymodis as pm
import re, os
from datetime import datetime, timedelta
from osgeo import ogr, gdal

def run_script(iface):
    
    #####################################################################################################################
    #####################################################################################################################
    ##PARAMETERS
    
    #Path to MRT software
    mrtPath = '/home/olivier/MRT'
    
    ##DIRECTORIES parameters
    #Working directory
    dst = '/media/olivier/olivier_ext/gedata_current/jde_coffee/MODIS/collection6' #'D:/gedata_current/jde_coffee/MODIS'
    #Folder to use for temporary files (should be empty)
    tempDir = 'Temp'
    #Destination folder for the download
    rawdataDir = 'raw_data'
    

    ##REGIONS parameters
    #Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    #Names of the regions (also folders names) 
    states = ["CER","CHA","CO","ES","MO","SDM","SP","ZM"]
    #Addresses of the shapefiles with the boundaries for each of the regions
    statesBoundFiles = ['/media/olivier/olivier_ext/gedata_current/jde_coffee/data/'+s+'/aoi/AOI_'+s+'.shp' 
                            for s in states] #Address of the boundary files !!!!SHOULD BE IN EPSG 4326 PROJECTION
    #Name of subfolder where to save the raw mosaic data (should be in the folders of the regions)
    statesRawFolder = 'raw_data'
    #Name of subfolder where to save the smoothed mosaic data (should be in the folders of the regions)
    statesSmoothFolder = 'smooth_data'
    
        
    ##DOWNLOAD parameters
    #Product to download
    product = 'MOD13Q1.006'
    #Username for the earthdata website
    user = "olivierpfrancois"
    #Password for the earthdata website
    pwd = "Michele1950"
    #Tiles to download
    tiles = ['h13v10','h13v11','h14v10','h14v11']
    #Start date for the product download (format YYYY-MM-DD)
    #    If None, will default to date of most recent MODIS file on disk if any, or stop the process
    startDownload = None
    #startDownload = '2017-04-10'
    #End date for the product download (format YYYY-MM-DD)
    #    If None, defaults to today
    endDownload = None
 
    
    ##MOSAIC parameters
    #Should the downloaded files be mosaiced for each of the regions?
    mosaic = 'no'
    #Starting date for the files to mosaic
    #    If None, will default to the files that have been just downloaded if any.
    startMosaic = None
    #startMosaic = '2017-04-10'
    #Ending date for the files to mosaic
    #    If None, defaults to today
    endMosaic = None
    
    
    ##SMOOTHING parameters
    smooth = 'no'
    
    
    #####################################################################################################################
    #####################################################################################################################
    ##FUNCTIONS
    def downloadMODIS(dstFolder, pwd, user, tiles, product, startDownload=None, endDownload=None):
        
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
        if not endMosaic:
            endMosaic = datetime.now()
            endMosaic = endMosaic.date()
            
        #Get the files on disk
        thereHdf = [f for f in os.listdir(os.path.join(root,srcFolder)) if f.endswith('.hdf')]
        #Dates of these files
        datesHdf = [re.search('A([0-9]{7})', f).group(1) for f in thereHdf]
        datesHdfD = [datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])) for d in datesHdf]
        
        #Get the dates to process
        dates = {d for d,D in zip(datesHdf,datesHdfD) if D >= startMosaic and D <= endMosaic}
        dates = sorted(list(dates))
        
        #Work by date
        for d in dates:
            #Transform into actual date
            dlong = datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:]))
            
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
            outName = re.sub('\.h[0-9]{2}v[0-9]{2}','',files[0])
            outName = re.sub('\.hdf','',outName)
            outName = re.sub('\.','_',outName)
            
            #Mosaic the files and transform to tif
            mos = pm.convertmodis_gdal.createMosaicGDAL(hdfnames=[os.path.join(root,srcFolder,f) for f in files], 
                        subset='1 0 0 0 0 0 0 0 0 0 0 0', outformat='GTiff') #outformat='HDF4Image'
            mos.run(output=os.path.join(root,tmpFolder,'temp.tif')) #'.hdf'
            #mos.write_mosaic_xml(os.path.join(root,tmpFolder,'outName'))
            #mos.write_vrt(os.path.join(root,tmpFolder,outName))
            
            #Change the projection to lat/long
            inRaster = '/'.join([root,tmpFolder,'temp.tif'])
            outRaster = '/'.join([root,tmpFolder,outName+'.tif'])
            cmd = 'gdalwarp -overwrite -t_srs EPSG:4326 -r near -of GTiff %s %s' % (inRaster, outRaster)
            os.system(cmd)

            '''
            src_ds = gdal.Open(os.path.join(root,tmpFolder,outName+'.hdf'))
            if not src_ds is None:
                layers = src_ds.GetSubDatasets()
                print(str(layers))
            conv = pm.convertmodis_gdal.convertModisGDAL(hdfname=os.path.join(root,tmpFolder,outName+'.hdf'), 
                            prefix=outName, subset='(1 0)', res=None, outformat='GTiff', 
                            epsg='4326', wkt=None, resampl='NEAREST_NEIGHBOR', vrt=False)
            conv.run()
            '''
            
            #Loop through the regions to clip and mask the mosaic
            for s, shp in zip(regions, regionsBoundaries):
                #Crop the mosaic and mask the raster
                inRaster = root+'/'+tmpFolder+'/'+outName+'.tif'
                outRaster = root+'/'+s+'/'+regionsOut+'/'+outName+'_ndvi.tif'
                cmd = 'gdalwarp -q -cutline %s -crop_to_cutline %s %s' % (shp, inRaster, outRaster)
                os.system(cmd)
            
            #Remove the .xml, .hdf and .vrt intermediary files
            xml = [f for f in os.listdir(os.path.join(root,tmpFolder)) if f.endswith('.xml') or 
                        f.endswith('.hdf') or f.endswith('.vrt') or f.endswith('.tif')]
            for f in xml:
                os.remove(os.path.join(root,tmpFolder,f))
                
    
    #####################################################################################################################
    #####################################################################################################################
    ##ACTIVE CODE
    
    newHdf = downloadMODIS(dstFolder=os.path.join(dst,rawdataDir), 
            pwd=pwd, user=user, tiles=tiles, product=product, startDownload=startDownload, 
            endDownload=endDownload)
    
    if mosaic == 'yes':
        
        thereHdf = [f for f in os.listdir(os.path.join(dst,rawdataDir)) if f.endswith('.hdf')]
        
        if not newHdf:
            if not thereHdf:
                return
            
            if not startMosaic:
                print('Provide at least a start date to mosaic files already on disk')
                return
        
        if startMosaic:
            startMosaic = datetime.strptime(startMosaic, '%Y-%m-%d').date()
        
        if endMosaic:
            endMosaic = datetime.strptime(endMosaic, '%Y-%m-%d').date()
                
        if newHdf:
            #Dates of these files
            datesNew = [re.search('A([0-9]{7})', f).group(1) for f in newHdf]
            datesNewD = [datetime.strptime('0101'+d[0:4], "%d%m%Y").date()+timedelta(days=int(d[4:])) for d in datesNew]
            if startMosaic:
                startMosaic = min(startMosaic, min(datesNewD))
            else:
                startMosaic = min(datesNewD)
         
        mosaicMODIS(root=dst, srcFolder=rawdataDir, tmpFolder=tempDir, 
                        regions=states, regionsOUt=statesRawFolder, regionsBoundaries=statesBoundFiles, tiles=tiles,
                        startMosaic=startMosaic, endMosaic=endMosaic)
    
    if smooth == 'yes':
        pass
        