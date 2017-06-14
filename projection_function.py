from osgeo import ogr, gdal
import osgeo.osr as osr

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



###########################################################################################################################################
##TESTING
        
##PARAMETERS
rName = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/ES/temp/01ES_201504191226247.TIF'

shpName1 = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/ES/temp/Espirito_Santo_geo.shp'
shpName2 = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/ES/temp/Espirito_Santo_utm.shp'
shpName3 = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/ES/temp/Espirito_Santo_utm_noproj.shp'


#Import the raster
raster = gdal.Open(rName)

#Import the shapefiles
driver = ogr.GetDriverByName('ESRI Shapefile')

shpgeo = driver.Open(shpName1, 1) # 0 means read-only. 1 means writeable.
lyrgeo = shpgeo.GetLayer(0)

shputm = driver.Open(shpName2, 1) # 0 means read-only. 1 means writeable.

shpno = driver.Open(shpName3, 1)

test = checkCompareProjections(src2=shpgeo, src1=raster, reproject=False)
#testl = test.GetLayer(0)
#print(test.GetProjection())
#print(testl.GetSpatialRef())
print(test)