# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" Your Description of the script goes here """

# Some commonly used imports

from datetime import datetime, timedelta
from osgeo import gdal, gdalconst, ogr
import os, re, sys, osr
import numpy as np

sys.path.append('/home/olivier/OlivierGithub/QGIS-scripts')
import projection_function as proj
import image_proc_functions as img

def run_script(iface):
    
    #Root directory for the modis data
    root = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6' #'E:/gedata_current/jde_coffee/MODIS/collection6'
    
    #Folder in root to use for temporary files (should be empty)
    tempDir = 'Temp'
    
    #Destination directory of the extracted modis information
    dst = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6'
    
    ##REGIONS parameters
    #Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    #Names of the regions (also folders names) 
    states = ['CER'] #["CER","CHA","CO","ES","MO","SDM","SP","ZM"]
    #Varieties in each case
    varieties = [['arabica']] #[['arabica'],['arabica'],['arabica'],['arabica','robusta'],['arabica'],['arabica'],['arabica'],['arabica','robusta']]
    #Name of subfolder where the smoothed images are located
    statesSmoothFolder = 'smooth_data'
    
    #Address of the shapefile with the polygons in which to extract
    shp = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/data/brazil_admin_boundaries/br_municipalities.shp'
    #Name of the field with the unique ids of the polygons
    idAttribute = 'CD_GEOCMU'
    
    #Start date for extraction
    startd = '2017-01-01'
    #End date for extraction
    endd = '2017-08-01'
    
    masks = [['masks/CER_densities_arabica_from_classifications_250m.tif']]
    '''
    [['masks/CER_densities_arabica_from_classifications_250m.tif'],
                    ['masks/CHA_densities_arabica_from_classifications_250m.tif'],
                    ['masks/CO_densities_arabica_from_classifications_250m.tif'],
                    ['masks/ES_densities_arabica_from_classifications_250m.tif','masks/ES_densities_robusta_from_classifications_250m.tif'],
                    ['masks/MO_densities_arabica_from_classifications_250m.tif'],
                    ['masks/SDM_densities_arabica_from_classifications_250m.tif'],
                    ['masks/SP_densities_arabica_from_classifications_250m.tif'],
                    ['masks/ZM_densities_arabica_from_classifications_250m.tif','masks/ZM_densities_robusta_from_classifications_250m.tif']]
    '''
    extractModis(root=root, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, outFolder=dst, masks=masks, 
                 startExtract=startd, endExtract=endd, shpName=shp, attr=idAttribute, tempDir=tempDir, epsgShp=4674)



def extractModis(root, regions, varieties, regionsIn, outFolder, masks, startExtract, endExtract, shpName, attr, tempDir, epsgShp=None):
    
    #Transform into date format
    if not endExtract:
        endExtract = datetime.now()
        endExtract = endExtract.date()
    else:
        endExtract = datetime.strptime(endExtract, '%Y-%m-%d').date()
    if not startExtract:
        startExtract = endExtract - timedelta(days=30)
    else:
        startExtract = datetime.strptime(startExtract, '%Y-%m-%d').date()
    
    
    #Open the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shpName, 0) # 0 means read-only. 1 means writeable.

    # Check if shapefile is correctly loaded.
    if dataSource is None: 
        print("ERROR","PROCESS ABORTED: Error (1) opening shapefile " + shpName)
        return False
    
    # Get layer
    inVecLayer = dataSource.GetLayer(0)
    
    #Get projection if possible
    inEpsg = proj.getProjection(dataSource)
    if not inEpsg and not epsgShp:
        print('Could not identify projection for the shapefile. Specify a value for epsgShp parameter')
        return False
    elif not inEpsg:
        inEpsg = epsgShp
    
    if not inEpsg == 4326:
        
        #Create the transform
        prj1 = osr.SpatialReference()
        prj1.ImportFromEPSG(inEpsg)
        prj2 = osr.SpatialReference()
        prj2.ImportFromEPSG(4326)
        coordTrans = osr.CoordinateTransformation(prj1,prj2)
        
        #create output datasource
        outdriver = ogr.GetDriverByName('ESRI Shapefile')
        srcOut = outdriver.CreateDataSource(os.path.join(root,tempDir,'temp.shp'))
        #Create a new layer
        lyrOut = srcOut.CreateLayer("layer 1", srs=prj2, geom_type=inVecLayer.GetGeomType())
        
        #Add fields
        lyr1Defn = inVecLayer.GetLayerDefn()
        for i in range(0, lyr1Defn.GetFieldCount()):
            fieldDefn = lyr1Defn.GetFieldDefn(i)
            lyrOut.CreateField(fieldDefn)
        
        # get the output layer's feature definition
        lyrOutDefn = lyrOut.GetLayerDefn()
        
        # loop through the input features
        inFeature = inVecLayer.GetNextFeature()
        while inFeature:
            # get the input geometry
            geom = inFeature.GetGeometryRef()
            # reproject the geometry
            geom.Transform(coordTrans)
            # create a new feature
            outFeature = ogr.Feature(lyrOutDefn)
            # set the geometry and attribute
            outFeature.SetGeometry(geom)
            for i in range(0, lyrOutDefn.GetFieldCount()):
                outFeature.SetField(lyrOutDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
            # add the feature to the shapefile
            lyrOut.CreateFeature(outFeature)
            # dereference the features and get the next input feature
            outFeature = None
            inFeature = inVecLayer.GetNextFeature()
        
        #Close the two shapefiles
        dataSource.Destroy()
        inVecLayer = None
        srcOut.Destroy()
        lyrOut = None
    
    #Loop through the regions to extract the modis data for each
    for r in range(len(regions)):
        print('Extracting region '+str(regions[r])+'...')
        
        #Import all the smooth modis images on disk
        onDisk = [f for f in os.listdir(os.path.join(root,regions[r],regionsIn)) if f.endswith('.tif')]
        
        #Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        #Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
        #Keep only the files and dates within the dates to process
        onDisk = [f for f,d in zip(onDisk,datesAll) if d >= startExtract and d <= endExtract]
        datesAll = [d for d in datesAll if d >= startExtract and d <= endExtract]
        #Sort the two list by date
        datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
        if not onDisk:
            print('no modis images to process in '+regions[r])
            continue
        
        for v in range(len(varieties[r])):
            #Get a base image as a reference for format
            baseImg = gdal.Open(os.path.join(root,regions[r],regionsIn,onDisk[0]))
            
            #Compute the sum for each polygon of the coffee
            coffeeWeights = img.raster_stats(regionsFileName=os.path.join(root,tempDir,'temp.shp'), rasterFileName=root+'/'+regions[r]+'/'+masks[r][v],  
                                             polIdFieldName=attr, statistics=['sum'], outTxt=None, addShp=False, addSuffix='', 
                                             numOutDecs=2, alltouch=True)
            
            if not coffeeWeights:
                print('error extracting coffee sums per polygon')
                continue
            
            coffeeWeights = coffeeWeights[0]
            
            #Import the mask
            weightsRaster = gdal.Open(root+'/'+regions[r]+'/'+masks[r][v])
            
            #nodataWeights = weightsRaster.GetRasterBand(1).GetNoDataValue()
            
            #Create an empty list to hold the averages by polygon for each date
            imgStat = []
            
            #Loop through the images
            for i, date in zip(onDisk,datesAll):
                print(date)
                
                #Import the image
                baseImg = gdal.Open(os.path.join(root,regions[r],regionsIn,i))
                
                #Remove existing temporary raster if any
                try:
                    os.remove(os.path.join(root,tempDir,'temp.tif'))
                except OSError:
                    pass
                
                #Create an empty copy for storing the temporary weighted image
                weighted = new_raster_from_base(baseImg, root+'/'+tempDir+'/temp.tif', 'GTiff', np.nan, gdal.GDT_Float32, bands=1)
                
                #Import image band as array
                imgBand = baseImg.GetRasterBand(1)
                imgArray = imgBand.ReadAsArray()
                
                #Recast the type to be sure
                imgArray = imgArray.astype(np.float32)
                #Replace the no data values by 0
                imgArray[np.logical_or(np.isnan(imgArray),np.logical_or(imgArray>1,imgArray<-1))] = 0.
                
                #Change the resolution of the raster to match the images if needed
                geoSmooth = weightsRaster.GetGeoTransform()
                geoBase = baseImg.GetGeoTransform()
                reproject = [1 for a,b in zip(geoSmooth,geoBase) if not a==b]
                if reproject:
                    weightsReproj = warp_raster(weightsRaster, baseImg, resampleOption='nearest', outputURI=None, outFormat='MEM')
                
                else:
                    weightsReproj = weightsRaster
                    
                #Import coffee weights band as array
                weightsBand = weightsReproj.GetRasterBand(1)
                weightsArray = weightsBand.ReadAsArray()
                weightsArray[weightsArray>1] = np.nan
                
                #Recast the type to be sure
                weightsArray = weightsArray.astype(np.float32)
                #Replace the no data values by 0
                weightsArray[np.logical_or(np.isnan(weightsArray),np.logical_or(weightsArray>1,weightsArray<-1))] = 0.
                
                #Estimate the weighted sum for each pixel
                weightedArray = np.multiply(weightsArray, imgArray)
                
                #Export to temporary raster
                weighted.GetRasterBand(1).WriteArray(weightedArray)
                
                #Close the temporary raster
                weighted.FlushCache()
                weighted = None
                
                stat = img.raster_stats(regionsFileName=os.path.join(root,tempDir,'temp.shp'), rasterFileName=root+'/'+tempDir+'/temp.tif',  
                                        polIdFieldName=attr, statistics=['sum'], outTxt=None, addShp=False, addSuffix='', 
                                        numOutDecs=2, alltouch=True)
                imgStat.append(stat[0])
                
            #Remove the temporary file
            try:
                os.remove(os.path.join(root,tempDir,'temp.tif'))
            except OSError:
                pass
            
            #Get the IDs that have some weighted ndvi value
            keys = set()
            for m in imgStat:
                keys.update(m.keys())
            
            #Combine the data from the coffee weights and the sum of the ndvi
            zz = 0
            for m in range(len(imgStat)):
                for k in imgStat[m].keys():
                    imgStat[m][k] = imgStat[m][k]['sum']/coffeeWeights[k]['sum']
                    if zz<3:
                        print(coffeeWeights[k]['sum'])
                        print(imgStat[m][k])
                        print(' ')
                        zz +=1
            
            #Extract the data for all the dates
        
        
        
                
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