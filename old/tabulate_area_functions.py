from osgeo import ogr, gdal
from csv import DictWriter
import os, processing, re
import numpy as np
import time
from qgis.core import *



def daan(a):
    
    print a


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
    outds.SetGeoTransform(geotransform) #data0.GetGeoTransform()  #get GeoTranform from existed 'data0'
    outds.SetProjection(projection) #data0.GetProjection()
    
    outds.FlushCache()
    outds=None
            
def GetExtentRaster(gt, cols, rows):
    ''' Return list of corner coordinates from a geotransform
        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    '''
    ext=[]
    xarr=[0,cols]
    yarr=[0,rows]
    
    #Goes through each of the four corners to compute the coordinates
    for px in xarr:
        for py in yarr:
            x=gt[0]+(px*gt[1])+(py*gt[2])
            y=gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
        yarr.reverse()
    return ext

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




def tabulate_area(rasterName, bands, shpName, index=None, tabulate=True, prefix='', statistics=None, outTxt=None, addShp=False, addSuffix='', alltouch=False):
    '''
    Computes the tabulate area for each polygon in the shapefile on each of the required bands
    
    INPUTS:
    rasterName: string
        Full address of the raster including extension
        
    bands: list of int
        list of band numbers to extract information from
    
    shpName: string
        Full address of the shapefile including extension
    
    index: [optional] string
        Name of the variable in the shapefile to use as labels for the polygons
    
    tabulate: [optional] logical
        If True the function will do a tabulate for each of the polygons
    
    prefix: [optional] string
        prefix to use to write in front of the raster values in the output file and dictionaries
    
    statistics: [optional] list of string
        The function will return for each polygon the statistics listed in this option. 
        The argument should be None if no statistics are wanted
        The statistics can be chosen among:
        'min','max','mean','count','sum','sd','median','unique'
        Where count is the number of pixels in the polygon, and 
        unique the number of unique values in the polygon
        THE STATISTICS ARE ONLY EXPORTED TO A TXT FILE.
    
    outTxt: [optional] string
        Full address of the output .txt including extension
    
    addShp: [optional] logical
        If True the summary statistics for each of the polygons will be added to the input shp file
    
    addSuffix: [optional] string
        Suffix to be used for the statistics names when exporting
    
    alltouch: [optional] logical
        Parameter value for the rasterization.
        If True, all pixels touched by the polygons will be considered as part of the polygons.
        If False, only those with the centroid inside the polygons will be considered.
        Defaults to False
        
    The function returns a list, one element for each band.
    Each element is a dictionary of dictionaries has the label of the polygon.
    The main dictionary has one element per polygon, where the polygon label is the key.
    The sub-dictionaries have one element per value found in the raster band, 
    where the value is the key and the count is the value.
    '''
    
    if not tabulate and not statistics:
        print("No tabulation or statistics specified.")
        return
    
    #Transform alltouch
    if alltouch:
        alltouch = 'TRUE'
    else:
        alltouch = 'FALSE'
    
    #Import the raster
    if isinstance(rasterName, str):
        try:
            ds = gdal.Open(rasterName)
        except RuntimeError, e:
            print 'Unable to open '+rasterName
            print e
            return
        
        if ds is None:
            print 'Error opening '+rasterName
            return
    #elif isinstance(rasterName, QgsRasterLayer):
    #    ds = rasterName
    else:
        print 'Input raster needs to be a file address' #or a QgsRasterLayer object
        return
    
    #Get the number of bands
    nBands = ds.RasterCount
    
    #Get the no data value
    no_data = ds.GetRasterBand(1).GetNoDataValue()
    
    #Get the data type for the array
    d_type = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
    if not d_type in ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]:
        print("The data type of the raster is not supported by the algorithm")
        return
    
    #Import the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shpName, 1) # 0 means read-only. 1 means writeable.
    # Check to see if shapefile is found.
    if dataSource is None:
        print 'Could not open %s' % (shpName)
        return
    else:
        layer = dataSource.GetLayer(0)
    
    #Add an id field for the rasterization
    layerNames = [layer.GetLayerDefn().GetFieldDefn(n).GetName() for n in range(layer.GetLayerDefn().GetFieldCount())]
    if not "IDRAST" in layerNames:
        idField = ogr.FieldDefn("IDRAST", ogr.OFTInteger)
        layer.CreateField(idField)
    idField = layer.GetLayerDefn().GetFieldIndex("IDRAST")
    
    #Get the index of the field to label values
    if index:
        #Get the index of the field with the polygon ids
        idx = layer.GetLayerDefn().GetFieldIndex(index)
    else:
        idx = idField
    
    #Prepare the labels for the output file
    labels = []
    i = 1
    for feat in layer:
        feat.SetField("IDRAST", i)
        layer.SetFeature(feat)
        i += 1
        labels.append(feat.GetFieldAsString(idx))
    #ids = range(1,i)
    '''
    #Get the extent of the shapefile layer
    shp_ext = layer.GetExtent()
    
    # Get Raster extent coordinates
    rast_ext = GetExtentRaster(ds.GetGeoTransform(), ds.RasterYSize, ds.RasterXSize)
    # Format to the same tuple structure as the vector layer
    rast_ext = (rast_ext[0][0],rast_ext[2][0],rast_ext[1][1],rast_ext[0][1])
    
    #Compare the two extent to see if raster inside the shapefile
    if (rast_ext[0] < shp_ext[0] or rast_ext[1] > shp_ext[1] or 
        rast_ext[2] < shp_ext[2] or rast_ext[3] > shp_ext[3]):
        
        d_types = ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]
        rtype = d_types.index(d_type)
        
        if isinstance(rasterName, str):
            clean = rasterName
            rasterName = QgsRasterLayer(rasterName, 'temp')
        
        QgsMapLayerRegistry.instance().addMapLayer(rasterName)
        
        #Crop the raster by the layer extent in case it is smaller
        args = {"INPUT":rasterName, "PROJWIN":", ".join(map(str, layer.GetExtent())), 
                    "RTYPE":rtype, "COMPRESS":0, "BIGTIFF":3, "OUTPUT":None}
        if no_data:
            args["NO_DATA"] = no_data
        
        ds = None
        clip = processing.runalg("gdalogr:cliprasterbyextent", args)
        ds = clip['OUTPUT'] #Get the address of the temporary file with the result of the algorithm
        ds = gdal.Open(ds) #Import again
        
        QgsMapLayerRegistry.instance().removeMapLayer(rasterName.id())
        try:
            os.remove(clean+'.aux.xml')
        except OSError:
            pass
        del rasterName
    '''
    
    #Rasterize the vector layer:
    #Prepare an empty raster to rasterize the shapefile
    #rVector = new_raster_from_base(ds, 'temp.tif', 'GTiff', -1, gdal.GDT_Int32)
    rVector = new_raster_from_base(ds, 'temp', 'MEM', -1, gdal.GDT_Int32)
    #Rasterize
    gdal.RasterizeLayer(rVector, [1], layer, options=['ALL_TOUCHED='+alltouch, 'ATTRIBUTE=IDRAST'])
    
    try:
        v_band = np.array(rVector.GetRasterBand(1).ReadAsArray())
        if len(v_band.shape) == 3:
            v_band = v_band[:,:,0] #Keep only the first two dimensions
        
    except ValueError:
        print("Raster file is too big for processing. Please crop the file and try again.")
        
        #Close the rasters
        ds = None
        rVector = None
        
        return
    
    #Close the raster
    rVector = None
    
    #Flatten the array
    v_band.ravel()
    
    #Get no data positions
    remove_nodataV = v_band!=-1
    
    #Remove the no data
    v_band = v_band[remove_nodataV]
    
    #Order the data
    order = np.argsort(v_band)
    v_band = v_band[order]
    
    #Prepare an empty list to hold the dictionaries for each band
    tabout = []
    
    #Prepare an empty list to hold the dictionaries for each band
    stats = []
    
    '''
    In each case the output is a list of dictionaries, with one element per band analyzed in the raster.
    In each dictionary the keys are the IDs of the polygons, and the values are dictionaries.
    For tabout the dictionary values are the values and counts, and for the stats dictionary
    they are the statistics names and values.
    '''
    
    if statistics:
        statistics = [s for s in statistics if s in ['min','max','mean','count','sum','sd','median','unique']]
        statistics.sort()
        
    for b in bands:
        #Convert the band to an array
        try:
            ds_band = np.array(ds.GetRasterBand(b).ReadAsArray())
            if len(ds_band.shape) == 3:
                ds_band = ds_band[:,:,0] #Keep only the first two dimensions
            
            #Flatten the array
            ds_band.ravel()
            
        except ValueError:
            print("Raster file is too big for processing. Please crop the file and try again.")
            
            #Close the rasters
            ds = None
            
            return
        
        #Remove no data values for the polygons
        ds_band = ds_band[remove_nodataV]
        
        #Order the data
        ds_band = ds_band[order]
        
        #Remove no data
        remove_nodataD = np.logical_and(np.logical_not(np.isnan(ds_band)), np.isfinite(ds_band))
        ds_band = ds_band[remove_nodataD]
        v_clean = v_band[remove_nodataD]
        
        #Get the indices of the polygon ids
        ids, pids = np.unique(v_clean, return_index=True)
        #pids = [np.searchsorted(v_clean, ind) for ind in ids]
        pids= pids.tolist()
        pids.append(-1)
        
        if tabulate:
            #Prepare a dictionary to hold the counts
            tabout.append({})
            
            for ind in range(len(ids)):
                
                #Get unique values and counts
                u, c = np.unique(ds_band[pids[ind]:pids[ind+1]], return_counts=True)
                
                #Store unique values and counts in tabulate dictionary
                tabout[-1][labels[ind]] = dict(zip([prefix+x for x in map(str,u)], c))
        
        if statistics:
            #Prepare a dictionary to hold the counts
            stats.append({})
            
            for ind in range(len(ids)):
                
                values = ds_band[pids[ind]:pids[ind+1]]
                values = values[values>=-1]
                res = {}
                
                #Get the statistics
                if 'min' in statistics:
                    res['min'] = np.nanmin(values)
                if 'max' in statistics:
                    res['max'] = np.nanmax(values)
                if 'mean' in statistics:
                    res['mean'] = np.nanmean(values, dtype=np.float64)
                if 'count' in statistics:
                    res['count'] = values.size
                if 'sum' in statistics:
                    res['sum'] = np.nansum(values, dtype=np.float64)
                if 'sd' in statistics:
                    res['sd'] = np.nanstd(values, dtype=np.float64)
                if 'median' in statistics:
                    res['median'] = np.nanmedian(values)
                if 'unique' in statistics:
                    u = np.unique(values)
                    res['unique'] = u.size
                
                #Store statistics in stats dictionary
                stats[-1][labels[ind]] = res
    
    if outTxt:
        if tabulate:
            for b, d in zip(bands, tabout):
                #Transform the result dictionary for the export
                out = []
                #Prepare a set to hold the unique values in the raster
                uniques = set()
                for k, v in d.items():
                    uniques.update(v.keys())
                    v['id'] = k
                    out.append(v)
                
                #Remove the no data value if any
                if no_data:
                    uniques.discard(n) 
                
                #Export to a tab delimited file
                if len(bands) > 1:
                    outNm = re.sub(".txt","_"+str(b)+".txt",outTxt)
                else:
                    outNm = outTxt
                
                with open(outNm, "w") as f:
                    dict_writer = DictWriter(f, ['id']+sorted(uniques, key = lambda x: float(re.sub('^'+prefix,'',x))), extrasaction='ignore', delimiter="\t", restval="0")
                    dict_writer.writeheader()
                    for p in out:
                        dict_writer.writerow(p)
        
        if statistics:
            for b, s in zip(bands, stats):
                #Transform the result dictionary for the export
                out = []
                #Prepare a set to hold the unique values in the raster
                for k, v in s.items():
                    v['id'] = k
                    out.append(v)
                
                #Export to a tab delimited file
                if len(bands) > 1:
                    outNm = re.sub(".txt","_STATS_"+str(b)+".txt",outTxt)
                    statNames = [w+addSuffix+str(b) for w in statistics]
                else:
                    outNm = re.sub(".txt","_STATS.txt",outTxt)
                    statNames = [w+addSuffix for w in statistics]
                    
                with open(outNm, "w") as f:
                    dict_writer = DictWriter(f, ['id']+statNames, extrasaction='ignore', delimiter="\t", restval="0")
                    dict_writer.writeheader()
                    for p in out:
                        dict_writer.writerow(p)
        
        
    if addShp and statistics:
        #Add the statistics to the shp as new variables
        
        #Get the field names for the shapepfile
        layerNames = [layer.GetLayerDefn().GetFieldDefn(n).GetName() for n in range(layer.GetLayerDefn().GetFieldCount())]
        
        for b, s in zip(bands, stats):
            #Prepare the field names
            if len(bands) > 1:
                statNames = [w+addSuffix+str(b) for w in statistics]
            else:
                statNames = [w+addSuffix for w in statistics]
            
            #Loop through the statistics to add the columns
            for nm in statNames:
                if not nm in layerNames:
                    idNew = ogr.FieldDefn(nm, ogr.OFTReal)
                    idNew.SetWidth(32)
                    idNew.SetPrecision(4) #added line to set precision
                    layer.CreateField(idNew)
            
            #Loop through the features to add the statistics values
            for feat in layer:
                for nm, w in zip(statNames, statistics):
                    #print nm+' '+w+' '+feat.GetFieldAsString(idx)+': '+str(s[feat.GetFieldAsString(idx)][w])
                    feat.SetField(nm, float(s[feat.GetFieldAsString(idx)][w]))
                layer.SetFeature(feat)
        
    #Remove rasterization id
    layer.DeleteField(idField)
    
    #Remove layer
    del layer
    del dataSource
    
    #Close the rasters
    ds = None
    
    '''
    #Clean the temporary files on disk
    try:
        os.remove(clip['OUTPUT'])
    except:
        pass
    '''
    #Return the list of dictionaries
    exports = {}
    if tabulate:
        exports['tabulate']=tabout
    if statistics:
        exports['statistics']=stats
    
    return exports
        
    
def ndvi_calc(rasterName, redBand, irBand, outTif=True, outName=None):
    '''
    Computes the ndvi from a raster and either exports it as a tif or simply returns it as gdal object
    
    INPUTS:
    rasterName: string
        Full address of the raster including extension
        
    redBand: int
        band number for the red band
    
    irBand: int
        band number for the IR band
    
    outTif: [optional] logical
        Defaults to True
        if True, the resulting ndvi raster will be exported to a .tif file
    
    outName: [optional] string
        Full address of the raster including extension.
        If no address is provided, the raster will be named adding a ndvi suffix to the input name

    '''    
        
        
    #Import the raster
    try:
        ds = gdal.Open(rasterName)
    except RuntimeError, e:
        print 'Unable to open '+rasterName
        print e
        return
    
    if ds is None:
        print 'Error opening '+rasterName
        return
    
    #Get the number of bands
    nBands = ds.RasterCount
    
    if nBands < max(redBand,irBand):
        print 'The band numbers provided exceed the number of bands in the input raster'
        return
    
    #Get the no data value
    no_data = ds.GetRasterBand(1).GetNoDataValue()
    
    #Get the data type for the array
    d_type = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
    if not d_type in ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]:
        print("The data type of the raster is not supported by the algorithm")
        return
    
    #Pass the two bands into arrays
    red = ds_band = np.array(ds.GetRasterBand(redBand).ReadAsArray())
    ir = ds_band = np.array(ds.GetRasterBand(irBand).ReadAsArray())
    red = red.astype(float)
    ir = ir.astype(float)
    
    #Replace the no data value 
    if no_data:
        red[red==float(no_data)] = np.nan
        ir[ir==float(no_data)] = np.nan
    
    #Compute the ndvi
    ndvi = (ir-red)/(ir+red)
    
    #Replace the no data value
    ndvi[np.isnan(ndvi)] = -99.
    
    if outTif:
        if not outName:
            outName = re.sub("(.tif|.TIF)", "_ndvi.tif", rasterName)
        
        #Check if raster already exists and delete it if is is the case
        try:
            os.remove(outName)
        except OSError:
            pass
        
        #Export the resulting raster
        rNdvi = new_raster_from_base(ds, outputURI=outName, format="GTiff", nodata=-99., datatype=gdal.GDT_Float32, bands=1)
        
        #Fill it with the values
        rNdvi.GetRasterBand(1).WriteArray(ndvi)
        
        rNdvi.FlushCache()
        rNdvi=None
        #Close raster
        ds = None
    else:
        #Export the resulting raster
        rNdvi = new_raster_from_base(ds, outputURI='ndvi', format="MEM", nodata=-99., datatype=gdal.GDT_Float32, bands=1)
        
        #Fill it with the values
        rNdvi.GetRasterBand(1).WriteArray(ndvi)
        
        #Close raster
        ds = None
        
        #Return the raster
        return rNdvi

def shpToRaster(shpName, field, resRaster):
    '''
    Rasterizes a vector layer at the wanted resolution
    
    shpName: string
        Full address of the shapefile including extension
    
    field: string
        Name of the field in the shapefile to be transferred in the raster
    
    resRaster: int/float
        Resolution for the output raster
    '''
    # Open the data source and read in the extent
    source_ds = ogr.Open(shpName)
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
    #Get the index of the field with the information for the raster
    definition = source_layer.GetLayerDefn()
    idx = definition.GetFieldIndex(field)
    #Get the type of the field
    fieldType = definition.GetFieldDefn(idx).GetType()
    fieldType = definition.GetFieldDefn(idx).GetFieldTypeName(fieldType)
    if fieldType == 'Integer':
        fieldType = gdal.GDT_Int16
    elif fieldType == 'Real':
        fieldType = gdal.GDT_Float32
    else:
        print('The field type is not covered by the algorithm')
    
    #Create the name of the output raster
    outName = re.sub(".shp",".tif",shpName)
    
    #Remove output file if it already exists
    try:
        os.remove(outName)
    except OSError:
        pass
    
    # Create the destination data source
    x_res = int((x_max - x_min) / resRaster)
    y_res = int((y_max - y_min) / resRaster)
    target_ds = gdal.GetDriverByName('GTiff').Create(outName, x_res, y_res, 1, fieldType, options=['COMPRESS=LZW'])
    target_ds.SetProjection(source_layer.GetSpatialRef().ExportToWkt())
    target_ds.SetGeoTransform((x_min, resRaster, 0, y_max, 0, -resRaster))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, options=['ALL_TOUCHED=FALSE', 'ATTRIBUTE='+field])
        
    target_ds = None
    source_ds = None

def rasterizeLayer(shpName, rasterName, indexVar, alltouch=False):
    '''
    Rasterizes a vector layer using a raster as an example
    
    shpName: string
        Full address of the shapefile including extension
    
    rasterName: string
        Full address of the raster including extension
    
    indexVar: string
        Name of the field in the shapefile to be used for the rasterization
    
    alltouch: [optional] logical
        Parameter value for the rasterization.
        If True, all pixels touched by the polygons will be considered as part of the polygons.
        If False, only those with the centroid inside the polygons will be considered.
        Defaults to False
    '''
    
    #Import the raster
    try:
        ds = gdal.Open(rasterName)
    except RuntimeError, e:
        print 'Unable to open '+rasterName
        print e
        return
    
    if ds is None:
        print 'Error opening '+rasterName
        return
    
    #Create an empty vector in memory
    rVector = new_raster_from_base(ds, re.sub(".shp",".tif",shpName), 'GTiff', -1, gdal.GDT_Int32)
    
    #Import the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shpName, 1) # 0 means read-only. 1 means writeable.
    # Check to see if shapefile is found.
    if dataSource is None:
        print 'Could not open %s' % (shpName)
        return
    else:
        layer = dataSource.GetLayer(0)
    
    #Rasterize the layer
    gdal.RasterizeLayer(rVector, [1], layer, options=['ALL_TOUCHED='+alltouch, 'ATTRIBUTE='+indexVar])