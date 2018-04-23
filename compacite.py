import numpy as np
from osgeo import gdal
import functools
import multiprocessing as mp
import gedata_tbox_raster as rt

from scipy.stats import _distn_infrastructure

'''
The algorithm for compactness takes a buffer around each missing value 
in the image and looks at the share of missing values in the buffer. 
If the missing values are very compact, they should have a large number 
of missing values surrounding them on average. If they are isolated 
they should have a small number every time.
'''

def compactnessWrapper(rastersIn, radius, nodata=None, parallel=False, nCores=None):
    if parallel and len(rastersIn)>1:
        if not nCores:
            nCores = mp.cpu_count()
        
        p = mp.Pool(nCores)
        
        pp = functools.partial(compactnessRaster, 
                               radius=radius, 
                               nodata=nodata)
        
        mpsComp = p.map(pp, rastersIn)
        
    else:
        mpsComp = []
        for r in rastersIn:
            mpsComp.append(compactnessRaster(rasterIn=r, radius=radius, 
                                     nodata=nodata))
        #mpsComp = [compactnessRaster(rasterIn=r, radius=radius, 
        #                             nodata=nodata) for r in rastersIn]
        
    return mpsComp
    

def compactnessRaster(rasterIn, radius, nodata=None):
    # Import that raster to get the positions of the missing values
    try:
        dst = gdal.Open(rasterIn)
    except RuntimeError: 
        return False
    
    #Fix the no data value
    if not nodata:
            nodata = dst.GetRasterBand(1).GetNoDataValue()
        
    if not nodata:
        return False
        
    if type(nodata) is int or nodata == np.nan:
        nodata = [nodata]
        
    # Import that raster to get the positions of the missing values
    arr = dst.GetRasterBand(1).ReadAsArray()
    
    #Get the no data value
    nodataOut = dst.GetRasterBand(1).GetNoDataValue()
    if not nodataOut:
        nodataOut = nodata[0]
        
    # Get the position of all the missing data
    mps = []
    
    for n in nodata:
        miss = np.where(arr == n)
        mps += map(list, zip(miss[0], miss[1]))
    
    # Loop through the missing data to replace in r
    mpsComp = [compactnessPixel(pixel=pixel, rasterDst=dst, 
                                nodataIn=nodata, radius=radius) 
                                for pixel in mps]
    
    #Compute the average share of missing pixels in surrounding buffers
    avg = sum(mpsComp)/float(len(mpsComp))
    
    return avg

def compactnessPixel(pixel, rasterDst, nodataIn, radius):
    
    #Get the bounds of the array to extract from the raster
    minX = max(0, pixel[0] - radius)
    maxX = pixel[0] + radius + 1
    minY = max(0, pixel[1] - radius)
    maxY = pixel[1] + radius + 1
    
    #Extract the array around the pixel
    arr = rt.readRasterBlock(src=rasterDst, colStart=minY,
                             rowStart=minX,
                             colBlockSize=maxY - minY,
                             rowBlockSize=maxX - minX,
                             band=1)
    
    #Get the position of the pixel in the extract
    pixX = min(pixel[0], radius)
    pixY = min(pixel[1], radius)
    
    #Create circular kernel to mask the values in the array
    k = createKernel(radius)
    
    #Cut the kernel to the correct dimensions
    k = k[(radius-pixX):,(radius-pixY):]
    k = k[:arr.shape[0],:arr.shape[1]]
        
    #Mask the array with the no data
    mask = np.isin(arr, nodataIn)
    
    #Mask the array with the kernel
    mask[k == 0] = False
    
    #Return the share of no data in the kernel
    return np.sum(mask)/float(np.sum(k))


# create circular kernel
def createKernel(radius):
    kernel = np.zeros((2*radius+1, 2*radius+1))
    y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    kernel[mask] = 1
    return kernel


if __name__ == '__main__':
    #Parallel option
    parallel = True
    nCores=3
    
    #Define the radius of pixels to consider for compactness
    rad = [3,5,7,10,15]
    
    #Define the value to consider as no data
    nodata = [-3000]
    
    root = ('/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/'+
                    'collection6/terra/Vietnam/LD/masked_missing/')
    
    #Test image to process
    rast = ['MOD13Q1_2018-02-02_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2018-01-17_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2018-01-01_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-12-19_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-12-03_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-11-17_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-11-01_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-10-16_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-09-30_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-09-14_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-08-29_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-08-13_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-07-28_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-07-12_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-06-26_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-06-10_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-05-25_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-05-09_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-04-23_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-04-07_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-03-22_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-03-06_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-02-18_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-02-02_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-01-17_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2017-01-01_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-12-18_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-12-02_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-11-16_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-10-31_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-10-15_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-09-29_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-09-13_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-08-28_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-08-12_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-07-27_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-07-11_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-06-25_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-06-09_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-05-24_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-05-08_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-04-22_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-04-06_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-03-21_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-03-05_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-02-18_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-02-02_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-01-17_h28v7_250m_16_days_NDVI.tif',
            'MOD13Q1_2016-01-01_h28v7_250m_16_days_NDVI.tif']
    
    import time
    for rr in rad:
        start = time.time()
        test = compactnessWrapper(rastersIn=[root+r for r in rast], radius=rr, nodata=-3000, 
                             parallel=parallel, nCores=nCores)
        print time.time() - start
        print test
    
    