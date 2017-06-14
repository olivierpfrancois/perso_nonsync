# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" Tool to split a shapefile of Ground Truth (GT) polygons into two rasters, 
    containing each half of the pixels of GT for each of the classes. 
    The script also produces a .txt file containing the labels of each 
    class for import into Grass.
    
    INPUTS:
    Each of the inputs is a list so that the script can work in batch mode.
    All lists should have the same number of elements. 
    The elements in the lists should be listed in the same order.
    
    ---satName (string): Full address of the raster image to be used as template
                    for the output rasters
    
    ---vtName (string): Full address of the shapefile with the GT polygons
    
    ---vtLabelsField (string): Name of the field in the GT shapefile with the 
                    labels to be used for the classifications
    
    ---outFolder (string): Address of the folder where to output the rasters and txt file
                    
    OUTPUTS:
    ---Two rasters (.tif) of the same name as the input shapefile, with suffix _train and _check
    ---One text file (.txt) of the same name as the input shapefile, with suffix _labels, containing 
        the correspondence between the raster values and the labels of the GT classes. 
        The values can also be found in the input shapefile in a new field class IDGRASS.
    
"""

# Imports
from osgeo import ogr, gdal
import numpy as np
import re, os
import osgeo.osr as osr

def run_script(iface):
    
    ########################################################################################################################
    ########################################################################################################################
    ##INPUT PARAMETERS
    #Separate the elements with a comma inside the lists if more than one element.
    #Elements in the different lists should be listed in the same order.
    
    satName = ['/media/olivier/olivier_ext/Temp/IMG_SPOT7_MS_20170417.TIF']
    
    vtName = ['/media/olivier/olivier_ext/Temp/VT_sheki_april_38N_V6.shp']
    
    vtLabelsField = ['C']
    
    outFolder = ['/media/olivier/olivier_ext/Temp']
    
    
    ########################################################################################################################
    ########################################################################################################################
    ##FUNCTIONS
    
    def arrayToRaster(array, outName, data_type, projection, geotransform):
        
        dims = array.shape
        if len(dims) == 2:
            y_pixels, x_pixels = dims
            bands = 1
        elif len(dims) == 3:
            y_pixels, x_pixels, bands = dims
            
        driver = gdal.GetDriverByName('GTiff')
        dataset = driver.Create(outName, x_pixels, y_pixels, bands, data_type)
        if bands == 1:
            dataset.GetRasterBand(1).WriteArray(array)
        else:
            for b in range(bands):
                dataset.GetRasterBand(b+1).WriteArray(array[:,:,b])
    
        #Add GeoTranform and Projection
        dataset.SetGeoTransform(geotransform)
        dataset.SetProjection(projection)
        
        dataset.FlushCache()
        dataset=None
    
    def new_raster_from_base(base, outputURI, rasterFormat, nodata, datatype, bands=None):
        ''' 
        Create an empty copy of a raster from an existing one
        
        base: gdal raster layer
            Name of the variable with the input raster to copy
        
        outputURI: string
            Address + name of the output raster (extension should agree with format, none for memory)
            
        rasterFormat: string
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
    
        driver = gdal.GetDriverByName(rasterFormat)
        
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
    
    def checkCompareProjections(src1, src2=None, reproject=False):
        '''
        If only src1, checks if projected and returns the epsg or None (reproject is ignored)
        
        If two sources:
            If reproject == True, returns the first source, reprojected to the projection of the second if necessary
                In case the first source is a shapefile, it creates a copy in memory of the dataSource with the new projection
                The input dataSource and shapefile on disk are therefore not modified.
                This is a different shapefile, so if modified in follow up code it should be exported explicitely
            If reproject == False, returns the epsg of the two sources as a tuple
        
        Returns None if one of the sources has no projection or if its epsg cannot be identified.
        
        src1, src2:
            if rasters, should be gdal Datasets (output from gdal.Open)
            if shapefiles, should be ogr DataSources (output from driver.Open)
        
        '''
        
        #Check the type of the inputs
        if type(src1) is gdal.Dataset:
            type1 = 'raster'
        elif type(src1) is ogr.DataSource:
            type1 = 'shp'
        else:
            raise ValueError('checkCompareProjections: First input should be gdal.Dataset or ogr.DataSource type')
        
        #Get the projection of the first source
        if type1 == 'raster':
            prj1 = src1.GetProjection() #Gets the projection in wkt
            if not prj1:
                raise ValueError('checkCompareProjections: First input has no projection')
            #Go from wkt to osr spatial reference
            prj1 = osr.SpatialReference(wkt=prj1)
        elif type1 == 'shp':
            #Get the layer from the data source
            lyr1 = src1.GetLayer(0)
            #Get the osr spatial reference
            prj1 = lyr1.GetSpatialRef()
            if not prj1:
                raise ValueError('checkCompareProjections: First input has no projection')
        
        #Extract the EPSG for the first input
        epsg1 = prj1.GetAttrValue('AUTHORITY', 1)
        if not epsg1: #If EPSG was not explicitly included in the projection info, try to guess if
            if prj1.AutoIdentifyEPSG() == 0: # success in identifying
                epsg1 = prj1.GetAuthorityCode(None) #Get the epsg
    
        
        if not src2 and not epsg1:
            raise ValueError('checkCompareProjections: First input EPSG could not be identified')
        elif not src2:
            #Export the epsg if it was found
            return epsg1
        
        
        #Get the type and projection of the second source
        if type(src2) is gdal.Dataset:
            type2 = 'raster'
        elif type(src2) is ogr.DataSource:
            type2 = 'shp'
        else:
            raise ValueError('checkCompareProjections: Second input should be gdal.Dataset or ogr.DataSource type')
        
        #Get the projection of the first source
        if type2 == 'raster':
            prj2 = src2.GetProjection() #Gets the projection in wkt
            if not prj2:
                raise ValueError('checkCompareProjections: Second input has no projection')
            #Go from wkt to osr spatial reference
            prj2 = osr.SpatialReference(wkt=prj2)
        elif type2 == 'shp':
            #Get the layer from the data source
            lyr2 = src2.GetLayer(0)
            #Get the osr spatial reference
            prj2 = lyr2.GetSpatialRef()
            if not prj2:
                raise ValueError('checkCompareProjections: Second input has no projection')
        
        #Extract the EPSG for the first input
        epsg2 = prj2.GetAttrValue('AUTHORITY', 1)
        if not epsg2: #If EPSG was not explicitly included in the projection info, try to guess if
            if prj2.AutoIdentifyEPSG() == 0: # success in identifying
                epsg2 = prj2.GetAuthorityCode(None) #Get the epsg
    
        if not epsg2:
            raise ValueError('checkCompareProjections: First input EPSG could not be identified')
        
        if not reproject:
            #Return the epsg codes for each of the inputs
            return (epsg1,epsg2)
        
        if epsg1 == epsg2:
            #Return the first input directly if there is no need to reproject
            return src1
        
        else:
            #Clean up the projections by getting them from the EPSG code
            prj1 = osr.SpatialReference()
            prj1.ImportFromEPSG(int(epsg1))
            prj2 = osr.SpatialReference()
            prj2.ImportFromEPSG(int(epsg2))
            
            if type1 == 'raster':
                dst_wkt = prj2.ExportToWkt() #Create wkt projection
                error_threshold = 0.125  #Error threshold --> use same value as in gdalwarp
                resampling = gdal.GRA_NearestNeighbour #Resampling method
                
                # Call AutoCreateWarpedVRT() to fetch default values for target raster dimensions and geotransform
                srcOut = gdal.AutoCreateWarpedVRT(src1,
                               None, # src_wkt : left to default value --> will use the one from source
                               dst_wkt,
                               resampling,
                               error_threshold)
            
            else:
                #Create the transform
                coordTrans = osr.CoordinateTransformation(prj1,prj2)
                
                #create output datasource in memory
                outdriver = ogr.GetDriverByName('MEMORY')
                srcOut = outdriver.CreateDataSource('memData')
                #Create a new layer
                lyrOut = srcOut.CreateLayer("layer 1", srs=prj2, geom_type=lyr1.GetGeomType())
                
                #Add fields
                lyr1Defn = lyr1.GetLayerDefn()
                for i in range(0, lyr1Defn.GetFieldCount()):
                    fieldDefn = lyr1Defn.GetFieldDefn(i)
                    lyrOut.CreateField(fieldDefn)
                
                # get the output layer's feature definition
                lyrOutDefn = lyrOut.GetLayerDefn()
                
                # loop through the input features
                inFeature = lyr1.GetNextFeature()
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
                    inFeature = lyr1.GetNextFeature()
            
            return srcOut
            
            
    
    
    ########################################################################################################################
    ########################################################################################################################
    ##ACTIVE CODE
    
    for k in range(len(satName)):
        
        #Import the satellite image that will serve as a template for the VT output
        if isinstance(satName[k], str):
            try:
                satRaster = gdal.Open(satName[k])
            except RuntimeError, e:
                print 'Unable to open '+satName[k]
                print e
                return
            
            if satRaster is None:
                print 'Error opening '+satName[k]
                return
        
        else:
            print 'Input raster needs to be a file address'
            return
        
        #Import the VT shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        vt = driver.Open(vtName[k], 1) # 0 means read-only. 1 means writeable.
        # Check to see if shapefile is found.
        if vt is None:
            print 'Could not open %s' % (vtName[k])
            return
        else:
            vtLayer = vt.GetLayer(0)
        
        #Check that the projection is the same for the inputs or reproject
        satRaster = checkCompareProjections(src1=satRaster, src2=vt, reproject=True)
        
        #Get the index of the field with the polygon Labels
        labelsIndex = vtLayer.GetLayerDefn().GetFieldIndex(vtLabelsField[k])
        
        #Get the unique values from the field with the VT classes names and create a numeric field for the rasterization
        labels = []
        for feat in vtLayer:
            #add all the values in the shapefile
            labels.append(feat.GetField(labelsIndex))
            #labels.append(feat.GetFieldAsString(labelsIndex))
        vtLayer.ResetReading()
        
        #Get the unique elements
        labels = set(labels)
        labels = list(labels)
        #Sort them 
        labels.sort()
        
        #Add a GRASS ID for each of the categories into a dictionary
        labels = {l: i+1 for l,i in zip(labels, range(len(labels)))}
        
        #Create a new field that will host the IDs for GRASS
        layerNames = [vtLayer.GetLayerDefn().GetFieldDefn(n).GetName() for n in range(vtLayer.GetLayerDefn().GetFieldCount())]
        if not "IDGRASS" in layerNames:
            idField = ogr.FieldDefn("IDGRASS", ogr.OFTInteger)
            vtLayer.CreateField(idField)
        
        #Add the IDs for GRASS
        for feat in vtLayer:
            feat.SetField("IDGRASS", labels[feat.GetField(labelsIndex)]) #Fill with the correct ID for that label
            vtLayer.SetFeature(feat) #Set the feature with this new value to save the change
        vtLayer.ResetReading()
        
        #Rasterize the VT layer based on the image
        #Prepare an empty raster to rasterize the shapefile
        #Link to the units http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4acd01d1afd29ffdb9a1cc6c09de61fd2b
        vtRaster = new_raster_from_base(satRaster, 'temp', 'MEM', 0, gdal.GDT_Byte, 1)
        #Rasterize
        gdal.RasterizeLayer(vtRaster, [1], vtLayer, options=['ALL_TOUCHED=TRUE', 'ATTRIBUTE=IDGRASS'])
        
        #Transform the VT raster into an array
        try:
            vtBand = np.array(vtRaster.GetRasterBand(1).ReadAsArray())
            if len(vtBand.shape) == 3:
                vtBand = vtBand[:,:,0] #Keep only the first two dimensions
            
        except ValueError:
            print("Raster file is too big for processing. Please crop the file and try again.")
            #Close the rasters
            satRaster = None
            vtRaster = None
        
        #Set half of the pixels to 0 for each class of the VT
        vtBandEven = vtBand.copy()
        vtBandOdd = vtBand.copy()
        for c in labels.itervalues():
            i = 0
            vtClass = np.where(vtBandEven == c)
            for index in zip(*vtClass):
                if i & 1:
                    vtBandOdd[index] = 0
                else:
                    vtBandEven[index] = 0
                i += 1
        
        #Export the resulting arrays
        outNm = outFolder[k] + '/' + os.path.basename(vtName[k])
        arrayToRaster(array=vtBandEven, outName=re.sub('.shp','_train.tif',outNm), 
                      data_type=gdal.GDT_Byte, projection=satRaster.GetProjection(), 
                      geotransform=satRaster.GetGeoTransform())
        arrayToRaster(array=vtBandOdd, outName=re.sub('.shp','_check.tif',outNm), 
                      data_type=gdal.GDT_Byte, projection=satRaster.GetProjection(), 
                      geotransform=satRaster.GetGeoTransform())
        
        #Export the labels into a txt file
        with open(re.sub('.shp','_labels.txt',outNm), 'w') as f:
            for k, v in labels.items():
                f.write('%s:%s\n' % (v, k))
