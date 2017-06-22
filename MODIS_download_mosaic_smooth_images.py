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
from osgeo import gdal
import numpy as np

def run_script(iface):
    
    #####################################################################################################################
    #####################################################################################################################
    ##PARAMETERS
    
    #Path to MRT software
    #mrtPath = '/home/olivier/MRT'
    
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
    #startDownload = '2017-05-26'
    #End date for the product download (format YYYY-MM-DD)
    #    If None, defaults to today
    endDownload = None
 
    
    ##MOSAIC parameters
    #Should the downloaded files be mosaiced for each of the regions?
    mosaic = 'yes'
    #Starting date for the files to mosaic
    #    If None, will default to the files that have been just downloaded if any.
    startMosaic = None
    #startMosaic = '2017-04-10'
    #Ending date for the files to mosaic
    #    If None, defaults to today
    endMosaic = None
    
    
    ##SMOOTHING parameters
    smooth = 'yes'
    #Starting date for the files to include as input in the smoothing process
    startSmooth = '2016-12-01'
    #Ending date for the files to include as input in the smoothing process
    #    If None, defaults to today
    endSmooth = None
    #Start and end dates for the files to save to the disk after smoothing
    startSaveSmooth = '2017-04-01' #None to save them all
    endSaveSmooth = None #None to save them up to the last date
    
    regWindow = 7
    avgWindow = 3
    
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
    
    
    def smoothMODIS(root, regions, regionsIn, regionsOut, startSmooth, endSmooth, startSaveSmooth=None, endSaveSmooth=None):
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
        startSaveSmooth (str): Starting date for the files to save. If None, will save all the processed files
        endSaveSmooth (str): Ending date for the files to save. If None, will save all the processed files
        '''
        
        startSmooth = datetime.strptime(startSmooth, '%Y-%m-%d').date()
        if not endSmooth:
            endSmooth = datetime.now()
            endSmooth = endSmooth.date()
        else:
            endSmooth = datetime.strptime(endSmooth, '%Y-%m-%d').date()
        
        #Loop through the regions to do the smoothing for each
        for r in regions:
            print('Processing region '+str(r)+'...')
            
            #Import all the raw modis images on disk
            onDisk = [f for f in os.listdir(os.path.join(root,r,regionsIn)) if f.endswith('.tif')]
            
            #Dates of these files
            datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
            #Transform into date format
            datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
            #Keep only the files within the dates to process
            onDisk = [f for f,d in zip(onDisk,datesAll) if d >= startSmooth and d <= endSmooth]
            onDisk.sort(key=dict(zip(onDisk, datesAll)).get) #Sort by date
            
            #Create a mask to know which of the dates to save
            if startSaveSmooth or endSaveSmooth:
                #Transform into date format
                if startSaveSmooth:
                    startSaveSmooth = datetime.strptime(startSaveSmooth, '%Y-%m-%d').date()
                if endSaveSmooth:
                    endSaveSmooth = datetime.strptime(endSaveSmooth, '%Y-%m-%d').date()
                #Sort the dates
                datesAll.sort()
                #Create a mask of which to save
                if startSaveSmooth and endSaveSmooth:
                    toSave = [d >= startSaveSmooth and d <= endSaveSmooth for d in datesAll]
                elif not startSaveSmooth:
                    toSave = [d <= endSaveSmooth for d in datesAll]
                elif not endSaveSmooth:
                    toSave = [d >= startSaveSmooth for d in datesAll]
            else:
                #Save them all
                toSave = [True for d in datesAll]
            
            #Import all the images to process in the right order
            toProcess = [gdal.Open(os.path.join(root,r,regionsIn,f)) for f in onDisk]
            
            #Get the no data value if any
            nodata = toProcess[1].GetRasterBand(1).GetNoDataValue()
            nodataCopy = nodata
            if not nodataCopy:
                nodataCopy = -1
            
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
                    print('   Processing block '+str(progress)+' of '+str(totBlocks))
                    progress += 1
                    
                    blocks = [read_raster_block(p, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize) for p in toProcess]
                    
                    #Bring the blocks together into one single array
                    blocks = np.dstack(blocks)
                    #Recast the type
                    blocks = blocks.astype(np.float32)
                    
                    #Do the smoothing for each pixel
                    blocks = smoothingSwets(blocks, regWindow, avgWindow, nodata)
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
                weights[0,0] = 0.5
                weights[len(pixel)-1,0] = 0.5
                for i in range(1, len(pixel)-1, 1):
                    if pixel[i-1,0] < pixel[i,0] and pixel[i,0] > pixel[i+1,0]:
                        weights[i,0] = 1.5
                    elif ((pixel[i-1,0] <= pixel[i,0] and pixel[i,0] <= pixel[i+1,0]) or (pixel[i-1,0] >= pixel[i,0] and pixel[i,0] >= pixel[i+1,0])):
                        weights[i,0] = 0.5
                    elif pixel[i-1,0] > pixel[i,0] and pixel[i,0] < pixel[i+1,0]:
                        weights[i,0] = 0.005
                
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
                    else:
                        res = np.mean(x[(int(np.floor(regWindow/2.))-t):(int(np.floor(regWindow/2.))+t+1)])    
                    smoothed.append(res)
                
                smoothed = np.asarray(smoothed)
                smoothed = np.reshape(smoothed, originShape)
                
                block[X,Y,:] = smoothed
        
        return block
    
    
    
    def read_raster_block(src, xStart, yStart, xBlockSize, yBlockSize, band=1):
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
        
        smoothMODIS(root=dst, regions=states, regionsIn=statesRawFolder, regionsOut=statesSmoothFolder, 
                    startSmooth=startSmooth, endSmooth=endSmooth, startSaveSmooth=startSaveSmooth, endSaveSmooth=endSaveSmooth)
        