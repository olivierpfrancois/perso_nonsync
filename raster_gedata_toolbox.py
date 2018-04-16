
from osgeo import gdal, gdalconst


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

def newRasterFromBase(base, outputURI, formatR, nodata, datatype, bands=None):
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
    datatype (gdal data type): (e.g. gdal.GDT_Int32) Data type for the raster
    bands [optional] (int): Number of bands for the output raster. 
        If not specified, will use the number of bands of the input raster
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

def readRasterBlock(src, colStart, rowStart, colBlockSize, rowBlockSize, band=1):
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
    
    # Get the wanted band
    band = src.GetRasterBand(band)
    
    # Get the size of the raster
    colSize = band.XSize
    rowSize = band.YSize
    
    rows = min(rowBlockSize, rowSize - rowStart)
    cols = min(colBlockSize, colSize - colStart)
            
    outArray = band.ReadAsArray(colStart, rowStart, cols, rows)
    
    return outArray

def warpRaster(src, dst, resampleOption='nearest', outputURI=None, outFormat='MEM'):
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
    
    # Raster to host the warped output
    if outFormat == 'MEM':
        rOut = newRasterFromBase(dst, 'temp', 'MEM', nodata, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    else:
        rOut = newRasterFromBase(dst, outputURI, outFormat, nodata, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    
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