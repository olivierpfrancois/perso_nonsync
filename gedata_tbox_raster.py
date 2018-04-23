
from osgeo import gdal, gdalconst
import numpy as np
import re

def createMask(src, values, masking, band=1, outFormat='MEM', outR=None):
    '''
    Fonction to create a mask from a raster dataset.
    The output raster will have values of 0 and 1, and 
    a no data value of 255
    
    src (gdal raster layer): Name of the variable with 
        the input raster to mask
    values (list): either a list of values to mask, 
        or a list of tuples that identify intervals 
        of values to mask/keep. 
        Tuples design closed intervals
    masking is one of two options, 'keep' or 'mask'.
        If 'keep', only the values provided will be 
        kept as 1. If 'mask', only these values will
        be passed to 0.
    bands (list): band to mask. 1 by default
    outFormat (str): format of output raster. If 'MEM', 
        will be created in memory and not saved to disk
        Other typical option is 'GTiff'
    outR (str): Full address of the output raster. Not 
        needed if outFormat is MEM
    '''
    
    #Check the input
    if not type(src) is gdal.Dataset:
        print('src should be a gdal dataset')
        return False
    
    # Raster to host the warped output
    if outFormat == 'MEM':
        rOut = newRasterFromBase(src, 'temp', 'MEM', 255, gdal.GDT_Byte, bands=1)
    else:
        rOut = newRasterFromBase(src, outR, outFormat, 255, gdal.GDT_Byte, bands=1)
    
    #Import the source raster as an array
    srcA = src.GetRasterBand(band).ReadAsArray()
    
    #Get the no data value if any
    nodata = src.GetRasterBand(1).GetNoDataValue()
    
    #Import the output raster as an array
    rOutA = rOut.GetRasterBand(band).ReadAsArray()
    
    #Switch all values to 0 or 1
    if masking == 'keep':
        if nodata:
            rOutA[srcA!=nodata] = 0
        else:
            rOutA = np.zeros(rOutA.shape, dtype=int)
    elif masking == 'mask':
        if nodata:
            rOutA[srcA!=nodata] = 1
        else:
            rOutA = np.ones(rOutA.shape, dtype=int)
    else:
        print('masking should be \'keep\' or \'mask\'')
        return False
    
    #Create the masked array
    if masking == 'keep':
        outVal = 1
    elif masking == 'mask':
        outVal = 0
    for v in values:
        if isinstance(v, int):
            rOutA[srcA==v] = outVal
        if isinstance(v, tuple):
            rOutA[srcA>=v[0] & srcA<=v[1]] = outVal
    
    rOut.GetRasterBand(band).WriteArray(rOutA)
    
    return rOut
    
    
def resampleRaster(src, downFactor, resampleOption='nearest', 
                     outputURI='', formatR='MEM', 
                     nodata=None, datatype=None):
    '''
    Fonction to resample a raster by a factor.
    The fonction creates an empty raster at the wanted resolution
    and then resamples the source raster into it using 
    gdal.ReprojectImage
    
    src (gdal Dataset): source raster to be warped
    downFactor (num): Factor by which to downsample the raster
        if >1, will downsample, if <1 will upsample
    resampleOption (string): One of 'nearest', 'bilinear', 'cubic', 
        'cubic spline', 'lanczos', 'average', or 'mode'. Method to use to 
        resample the pixels of the source raster
    outputURI (string, optional): Full address and name of the output raster. 
        If outFormat is 'MEM',this argument is ignored and the function simply 
        produces a raster in memory. 
        The extension for the output file should match the outFormat.
    formatR (string, optional): Format to use for the output raster from 
        the function. 
        Use 'GTiff' for a .tif output. Default creates a raster in memory.
    nodata (num): No data value (type should agree with raster type)
        If none is provided, it will use the input no data if any
    datatype (gdal const): data type for the output raster 
        (e.g. gdal.GDT_Float32).
        If None, will use th data type of the input raster
    '''
    
    #Get the raster information
    cols = src.RasterXSize
    rows = src.RasterYSize
    projection = src.GetProjection()
    geoT = src.GetGeoTransform()
    bands = src.RasterCount
    if not nodata:
        nodata = src.GetRasterBand(1).GetNoDataValue()
    if not nodata:
        print('No no data information in input raster \
                and none provided')
        return False
    if not datatype:
        datatype = src.GetRasterBand(1).DataType
    
    #Create a raster at lower resolution
    driver = gdal.GetDriverByName(formatR)
    if not outputURI and not formatR == 'MEM':
        print('Need to provide an output name for that format')
        return False
        
    if formatR == "GTiff":
        rOut = driver.Create(outputURI, 
                             int(cols/downFactor), 
                             int(rows/downFactor), 
                             bands, 
                             datatype, 
                             options=['COMPRESS=LZW', 'BIGTIFF=IF_NEEDED'])
    else:
        rOut = driver.Create(outputURI, 
                             int(cols/downFactor), 
                             int(rows/downFactor), 
                             bands, 
                             datatype)
    rOut.SetProjection(projection)
    #Create new geotransform with the correct resolution
    outGeoT = (geoT[0], geoT[1]*downFactor, geoT[2], geoT[3], \
            geoT[4], geoT[5]*downFactor)
    rOut.SetGeoTransform(outGeoT)
    
    #Fill the bands with no data
    for i in range(bands):
        rOut.GetRasterBand(i + 1).SetNoDataValue(nodata)
        rOut.GetRasterBand(i + 1).Fill(nodata)
    
    # Define resampling options
    resampleOptions = {'nearest': gdalconst.GRA_NearestNeighbour, 'bilinear':gdalconst.GRA_Bilinear,
                   'cubic':gdalconst.GRA_Cubic, 'cubic spline':gdalconst.GRA_CubicSpline,
                   'lanczos':gdalconst.GRA_Lanczos, 'average':gdalconst.GRA_Average,
                   'mode':gdalconst.GRA_Mode} 
    
    if not resampleOption in resampleOptions.keys():
        print('Resampling option provided is not allowed')
        return False
    
    # Warp: the parameters are source raster, destination raster, 
    #        source projection, destination projection, 
    #        resampling option 
    gdal.ReprojectImage(src, rOut, projection, projection, 
                        resampleOptions[resampleOption])
    
    return rOut
    
    
def getGDALTypeFromNumber(nb):
    '''
    Takes the number representing the gdal data type and returns the actual
    data type as a string.
    '''
    
    if nb == 0:
        out = 'GDT_Unknown'
    elif nb == 1:
        out = 'GDT_Byte'
    elif nb == 2:
        out = 'GDT_UInt16'
    elif nb == 3:
        out = 'GDT_Int16'
    elif nb == 4:
        out = 'GDT_UInt32'
    elif nb == 5:
        out = 'GDT_Int32'
    elif nb == 6:
        out = 'GDT_Float32'
    elif nb == 7:
        out = 'GDT_Float64'
    elif nb == 8:
        out = 'GDT_CInt16'
    elif nb == 9:
        out = 'GDT_CInt32'
    elif nb == 10:
        out = 'GDT_CFloat32'
    elif nb == 11:
        out = 'GDT_CFloat64'
    elif nb == 12:
        out = 'GDT_TypeCount'
    
    return(out)

def newRasterFromBase(base, outputURI, formatR, nodata=None, 
                      datatype=None, bands=None):
    ''' 
    Create an empty copy of a raster from an existing one.
    The function returns a gdal raster variable filled with the nodata 
    value.
    
    base (gdal raster layer): Name of the variable with the input raster 
        to copy
    outputURI (str): Address + name of the output raster (extension 
        should agree with format, none for memory)
    formatR (str): Format for the dataset (e.g. "GTiff", "MEM")
    nodata (num): No data value (type should agree with raster type)
        If none is provided, it will use the input no data if any
    datatype (gdal data type): (e.g. gdal.GDT_Int32) Data type for the raster.
        If None is provided, it will use the input raster data type
    bands [optional] (int): Number of bands for the output raster. 
        If not specified, will use the number of bands of the input raster
    '''
    
    cols = base.RasterXSize
    rows = base.RasterYSize
    projection = base.GetProjection()
    geotransform = base.GetGeoTransform()
    if not bands:
        bands = base.RasterCount
    if not nodata:
        nodata = base.GetRasterBand(1).GetNoDataValue()
    if not nodata:
        print('No no data information in input raster \
                and none provided')
        return False
    
    if not datatype:
        datatype = base.GetRasterBand(1).DataType
        
    driver = gdal.GetDriverByName(formatR)
    
    if formatR == "GTiff":
        new_raster = driver.Create(str(outputURI), cols, rows, 
                                   bands, datatype, 
                                   options=['COMPRESS=LZW', 'BIGTIFF=IF_NEEDED'])
    else:
        new_raster = driver.Create(str(outputURI), cols, rows, 
                                   bands, datatype)
    new_raster.SetProjection(projection)
    new_raster.SetGeoTransform(geotransform)

    for i in range(bands):
        new_raster.GetRasterBand(i + 1).SetNoDataValue(nodata)
        new_raster.GetRasterBand(i + 1).Fill(nodata)

    return new_raster


def readRasterBlock(src, colStart, rowStart, colBlockSize, 
                    rowBlockSize, band=1):
    '''
    Function to read a block of data from a gdal Dataset. 
    It will return the block as a numpy array
    
    src (gdal Dataset): raster to extract the block from
    xStart (int): X pixel value from which to start extraction (upper left corner)
    yStart (int): Y pixel value from which to start extraction (upper left corner)
    xBlockSize (int): X size of the block to extract. If the block is larger than the number 
            of pixels from xStart, will extract up to the raster limit
    yBlockSize (int): Y size of the block to extract. If the block is larger than the number 
            of pixels from xStart, will extract up to the raster limit
    band (int): band from which to extract the block. 1 by default.
    '''
    
    # Get the wanted band
    band = src.GetRasterBand(band)
    
    # Get the size of the raster
    colSize = band.XSize
    rowSize = band.YSize
    
    rows = min(rowBlockSize, rowSize - rowStart)
    cols = min(colBlockSize, colSize - colStart)
            
    outArray = band.ReadAsArray(colStart, rowStart, cols, rows)
    
    return outArray

def warpRaster(src, dst, resampleOption='nearest', outDataType=None, outputURI=None, outFormat='MEM'):
    '''
    Function : Warp a source raster to the resolution, extent and projection of a destination raster.
           
    The function returns the resulting raster. If outFormat is different from 'MEM', the 
    raster is also saved to disk using the information provided in outputURI.
    
    
    src (gdal Dataset): source raster to be warped
    dst (gdal Dataset): destination raster that will provide the resolution, 
        extent and projection
    resampleOption (string): One of 'nearest', 'bilinear', 'cubic', 
        'cubic spline', 'lanczos', 'average', or 'mode'. Method to use to 
        resample the pixels of the source raster
    outDataType (gdal const): data type for the output raster 
        (e.g. gdal.GDT_Float32).
        If None, will use th data type of the input raster
    outputURI (string, optional): Full address and name of the output raster. 
        If outFormat is 'MEM',this argument is ignored and the function simply 
        produces a raster in memory. 
        The extension for the output file should match the outFormat.
            
    outFormat (string, optional): Format to use for the output raster from 
        the function. 
        Use 'GTiff' for a .tif output. Default creates a raster in memory.
    '''
    
    if not type(src) is gdal.Dataset:
        return False
    
    if not type(dst) is gdal.Dataset:
        return False
    
    # Define resampling options
    resampleOptions = {'nearest': gdalconst.GRA_NearestNeighbour, 'bilinear':gdalconst.GRA_Bilinear,
                   'cubic':gdalconst.GRA_Cubic, 'cubic spline':gdalconst.GRA_CubicSpline,
                   'lanczos':gdalconst.GRA_Lanczos, 'average':gdalconst.GRA_Average,
                   'mode':gdalconst.GRA_Mode} 
    
    if not resampleOption in resampleOptions.keys():
        return False
    
    nodata = src.GetRasterBand(1).GetNoDataValue()
    if not nodata:
        nodata = 0
    
    if not outDataType:
        outDataType = src.GetRasterBand(1).DataType
    
    # Raster to host the warped output
    if outFormat == 'MEM':
        rOut = newRasterFromBase(dst, 'temp', 'MEM', nodata, outDataType, bands=src.RasterCount)
    else:
        rOut = newRasterFromBase(dst, outputURI, outFormat, nodata, outDataType, bands=src.RasterCount)
    
    # Warp: the parameters are source raster, destination raster, source projection, destination projection, resampling option 
    gdal.ReprojectImage(src, rOut, src.GetProjection(), rOut.GetProjection(), resampleOptions[resampleOption])
    
    return rOut

def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    Returns values based on origin and resolution of the raster. 
    The actual values can fall outside of the raster and should be checked.
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    col = int((x - ulX) / xDist)
    row = int((ulY - y) / xDist)
    
    return (col, row)