import numpy as np
from osgeo import ogr, gdal


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



inr = gdal.Open('/media/olivier/olivier_ext/gedata_current/jde_coffee/MODIS/test/ES/raw_data/MOD13Q1_2016-12-02_h13-14v10-11_250m_16_days_NDVI.tif')

datatype = inr.GetRasterBand(1).DataType

cols = inr.RasterXSize
rows = inr.RasterYSize
projection = inr.GetProjection()
geotransform = inr.GetGeoTransform()
bands = inr.RasterCount
driver = gdal.GetDriverByName('GTiff')
outr = driver.Create(str('/media/olivier/olivier_ext/gedata_current/jde_coffee/MODIS/test/ES/raw_data/test.tif'), cols, rows, bands, datatype, options=['COMPRESS=LZW'])
outr.SetProjection(projection)
outr.SetGeoTransform(geotransform)

#outr = new_raster_from_base(inr, '/media/olivier/olivier_ext/gedata_current/jde_coffee/MODIS/test/ES/raw_data/test.tif', 'GTiff', -1, datatype)

band = inr.GetRasterBand(1).ReadAsArray()
outrBand = outr.GetRasterBand(1)
outrBand.writeArray(band)

block = read_raster_block(inr, 256, 256, 256, 256)
outr.GetRasterBand(1).writeArray(np.reshape(block, (256,256)), 256, 256)



