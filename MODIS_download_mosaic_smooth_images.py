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
import re, os
from datetime import datetime, timedelta
from osgeo import ogr, gdal
import numpy as np

def run_script(iface):
    
    #####################################################################################################################
    #####################################################################################################################
    ##PARAMETERS
    
    #Path to MRT software
    mrtPath = '/home/olivier/MRT'
    
    ##DIRECTORIES parameters
    #Working directory
    dst = '/media/olivier/olivier_ext/gedata_current/jde_coffee/MODIS/test' #'D:/gedata_current/jde_coffee/MODIS'
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
    #startDownload = None
    startDownload = '2017-05-26'
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
    smooth = 'yes'
    startSmooth = '2016-12-01'
    endSmooth = '2017-05-30'
    
    
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
    
    
    def smoothMODIS(root, regions, regionsIn, regionsOut, startSmooth, endSmooth):
        #Loop through the regions to do the smoothing for each
        for r in regions:
            
            #Import all the raw modis images on disk
            onDisk = [f for f in os.listdir(os.path.join(root,r,regionsIn)) if f.endswith('.tif')]
            #Dates of these files
            datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
            #Transform into date format
            datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
            #Keep only the files within the dates to process
            onDisk = {f for f,d in zip(onDisk,datesAll) if d >= startSmooth and d <= endSmooth}
            
            #Import all the images to process
            toProcess = [gdal.Open(os.path.join(root,r,regionsIn,f)) for f in onDisk]
            
            #Get the data type
            typeData = toProcess[1].GetRasterBand(1).DataType
            
            #Create empty copies of these files to use for the smoothed data
            processed = [new_raster_from_base(p, root+'/'+r+'/'+regionsOut+'/smooth_'+f, 'GTiff', -1, typeData) for p,f in zip(toProcess,onDisk)]
            
            #Get the size of the rasters to identify the limits of the blocks to loop through
            band = toProcess[1].GetRasterBand(1)
            #Get the size of the raster
            xsize = band.XSize
            ysize = band.YSize
            #Get the number of blocks in x and y directions
            xBlocks = int(round(xsize/256)) + 1
            yBlocks = int(round(ysize/256)) + 1
            
            for xStep in range(xBlocks):
                for yStep in range(yBlocks):
                    
                    blocks = [read_raster_block(p, xStep*256, yStep*256, 256, 256) for p in toProcess]
                    
                    #Bring the blocks together into one single array
                    blocks = np.dstack(blocks)
                    
                    #Do the smoothing for each pixel
                    
                    
                    #Change the values in the output raster
                    blocks = np.dsplit(blocks, len(toProcess))
                    for p, b in zip(processed, blocks):
                        p.GetRasterBand(1).writeArray(np.reshape(b, (256,256)), xStep*256, yStep*256)
                    #processed = [p.GetRasterBand(1).writeArray(np.reshape(b, (256,256)), xStep*256, yStep*256) for p,b in zip(processed, blocks)]
            
            #Close the rasters
            for p in processed:
                p.FlushCache()
                p = None
            
            for p in toProcess:
                p=None
    
    def smoothingSwets(pixel, regWindow, avgWindow):
        pass
        '''
        #Check if vector is NA
        if (length(which(!is.na(raw)))==0) {
          smooth <- as.numeric(rep(NA,length(raw)))
        } else {
          #Interpolate missing values 
          #THIS SHOULD BE DONE BEFOREHAND FOR THE ENTIRE RASTER USING TEMPORAL AND SPATIAL INFERENCE
          z <- which(is.na(raw))
          nz <- length(z)
          nx <- length(raw)
          if (nz > 0 & nz < nx) { 
            raw[z] <- spline(x=1:nx, y=raw, xout=z, method="natural")$y
          }
          #Create vector of weights for regressions
          weights <- 0.5 #Consider the first point to be a local sloping point
          for (i in 2:(length(raw)-1)) { #Local peak = 1.5, local valley = 0.005, local sloping = 0.5
            if (raw[i-1] < raw[i] & raw[i] > raw[i+1]) {
              weights <- c(weights,1.5)
            } else if ((raw[i-1] <= raw[i] & raw[i] <= raw[i+1]) | (raw[i-1] >= raw[i] & raw[i] >= raw[i+1])) {
              weights <- c(weights,0.5)
            } else if (raw[i-1] > raw[i] & raw[i] < raw[i+1]) {
              weights <- c(weights,0.005)
            }
          }
          weights <- c(weights,0.5) #Consider the last point to be a sloping point as well
          
          #Create a matrix with the ndvi data for this pixel and for the weights
          #For the data:
          #Each column will be the same
          ndvi.raw <- matrix(data=raw, nrow=length(raw), ncol=(length(raw)+reg.window))
          #Set to NA all data for each column besides the data used for the regressions
          ndvi.raw[lower.tri(ndvi.raw)] <- NA
          ndvi.raw[,reg.window:ncol(ndvi.raw)][upper.tri(ndvi.raw[,reg.window:ncol(ndvi.raw)])] <- NA
          #Remove the first two and last 3 columns, since they don't have enough points for a regression
          ndvi.raw <- ndvi.raw[,3:(ncol(ndvi.raw)-3)]
          #For the weights:
          weights <- matrix(data=weights, nrow=nrow(ndvi.raw), ncol=ncol(ndvi.raw))
          weights[which(is.na(ndvi.raw))] <- NA
          
          #Create regression results matrix
          ndvi.reg <- ndvi.raw
          #Estimate the regressions for each column
          for (i in 1:ncol(ndvi.raw)) {
            #Prepare data
            y <- ndvi.raw[,i][which(!is.na(ndvi.raw[,i]))]
            w <- weights[,i][which(!is.na(ndvi.raw[,i]))]
            x <- c(1:length(y))
            #Estimate potential outliers
            out <- abs(y-mean(y))/sd(y) #Will be an outlier if greater than 3 (only extreme outliers should be picked up)
            #Remove outliers before regression
            y.out <- y[which(out<3)]
            w.out <- w[which(out<3)]
            x.out <- x[which(out<3)]
            #Compute parameters of regression
            b <- (sum(w.out)*sum(w.out*x.out*y.out)-sum(w.out*y.out)*sum(w.out*x.out))/
              (sum(w.out)*sum(w.out*x.out^2)-sum(w.out*x.out)^2)
            a <- (sum(w.out*y.out)-b*sum(w.out*x.out))/sum(w.out)
            #Compute predicted values from regression
            ndvi.reg[,i][which(!is.na(ndvi.raw[,i]))] <- a + b*x
          }
          #Combination of the results from the regression for each point
          #Now we combine for each point the results from the moving regression windows
          #We take into account the results from avg.window regression windows, centered around the point, unless we are at the edges
          smooth <- c() #Will hold the smoothed results
          t <- floor(avg.window/2) #number of predicted values from regressions to take into account around the center
          for (i in 1:nrow(ndvi.reg)) {
            x <- ndvi.reg[i,][which(!is.na(ndvi.reg[i,]))]
            if ((i < ceiling(reg.window/2)) & (length(x) < reg.window)) {
              res <- 
                mean(x[max(1,ceiling(reg.window/2)-(reg.window-length(x))-t):max(2,ceiling(reg.window/2)-(reg.window-length(x))+t)])
            } else {
              res <- mean(x[(ceiling(reg.window/2)-t):(ceiling(reg.window/2)+t)])
            }
            smooth <- c(smooth,res)
          }
        }
        #ndvi <- array(data=ndvi,dim=c(1,1,length(x)))
        return(smooth)
        '''
    
    def read_raster_block(src, xStart, yStart, xBlockSize, yBlockSize, band=1):
        
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
            
    def new_raster_from_base(base, outputURI, format, nodata, datatype, bands=None):
        ''' 
        Create an empty copy of a raster from an existing one
        
        base: gdal raster layer
            Name of the variable with the input raster to copy
        
        outputURI: string
            Address + name of the output raster (extension should agree with format, none for memory)
            
        format: string
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
    
        driver = gdal.GetDriverByName(format)
        
        if format == "GTiff":
            new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype, options=['COMPRESS=LZW'])
        else:
            new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype)
        new_raster.SetProjection(projection)
        new_raster.SetGeoTransform(geotransform)
    
        for i in range(bands):
            new_raster.GetRasterBand(i + 1).SetNoDataValue(nodata)
            new_raster.GetRasterBand(i + 1).Fill(nodata)
    
        return new_raster            
    
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
        if not startSmooth or not endSmooth:
            print('The start and end dates are necessary for the smoothing')
            return
        
        startSmooth = datetime.strptime(startSmooth, '%Y-%m-%d').date()
        endSmooth = datetime.strptime(endSmooth, '%Y-%m-%d').date()
        
        smoothMODIS(root=dst, regions=states, regionsIn=statesRawFolder, regionsOut=statesSmoothFolder, startSmooth=startSmooth, endSmooth=endSmooth)
        