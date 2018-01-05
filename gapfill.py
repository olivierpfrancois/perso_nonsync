# Implements the gapfill algorithm as in the gapfill package in R

import numpy as np
from datetime import datetime
from osgeo import gdal
import MODIS_gedata_toolbox as md
from statsmodels.distributions.empirical_distribution import ECDF
import scipy.stats as ss
import pandas as pd
import statsmodels.formula.api as smf
from xml.dom import NoDataAllowedErr
import MODIS_gedata_toolbox as md
import os, re


def array2TwoDim(a):
    # Transforms the output of array2Matrix into a numpy array in 2 dim
    
    s = a.shape
    
    return(a.reshape(s[0] * s[1], s[2]))


def arrayAroundRasters(rasters, seasons, years, mp, size, nodata=None):
    '''
    Takes a subset of size size around a pixel located in position mp from
    a list of rasters
    
    rasters (list): list of rasters already opened in gdal
    seasons (list): Same size as rasters, gives for each the season
    years (list): Same size as rasters, gives for each the year
    mp (list): list giving the position of the pixel around which to extract.
        The first two elements are x and y, the third is the season and the 
        fourth is the year.
    size (list): size of the array to extract around the center pixel.
        The first two elements are x and y, the third is the season and the 
        fourth is the year.
    nodata (list): list of values to be considered as no data, or None.
        If list, all no data values will be changed for the first one. 
    
    returns a numpy array with the size dimension
    '''
    
    # Get the max and min positions inside the rasters
    minX = max(0, mp[0] - size[0])
    maxX = mp[0] + size[0] + 1
    minY = max(0, mp[1] - size[1])
    maxY = mp[1] + size[1] + 1
    
    # Get the unique season values, sorted
    S = set(seasons)
    S = list(S)
    S.sort()
    
    # Get the max and min season positions
    minSeason = max(0, S.index(mp[2]) - size[2])
    maxSeason = min(len(S) - 1, S.index(mp[2]) + size[2] + 1)
    
    # Get the unique year values, sorted
    Y = set(years)
    Y = list(Y)
    Y.sort()
    
    # Get the max and min year positions
    minYear = max(0, Y.index(mp[3]) - size[3])
    maxYear = min(len(Y) - 1, Y.index(mp[3]) + size[3] + 1)
    
    # Subset the list and get the position of the raster of interest
    subset = []
    i = 0
    for r, s, y in zip(rasters, seasons, years):
        # Add the raster
        if (s in range(S[minSeason], S[maxSeason] + 1) and 
            y in range(Y[minYear], Y[maxYear] + 1)):
            subset.append(r)
            # Update the position of the raster of interest in the subset
            if y < mp[3]:
                i += 1
            elif y == mp[3]:
                if s < mp[2]:
                    i += 1
    
    block = [md.readRasterBlock(src=r,
                                colStart=minY,
                                rowStart=minX,
                                colBlockSize=maxY - minY,
                                rowBlockSize=maxX - minX,
                                band=1) for r in subset]
    
    block = np.dstack(block)
    
    # Change the no data values
    if nodata and len(nodata) > 1:
        for n in nodata[1:]:
            block[block == n] = nodata[0]
    
    # Get the position of the pixel inside the extract
    pos = [
        min(mp[0], size[0]),
        min(mp[1], size[1]),
        i]
    
    return [block, pos]


def arrayAroundArray (a, mp, size):
    
    # Get the max and min positions inside the array
    minPos = [max(0, m - s) for m, s in zip(mp, size)]
    maxPos = [min(dim, m + s + 1)for dim, m, s in zip(a.shape, mp, size)]
    
    # Transform into ranges
    n = [range(m, M) for m, M in zip(minPos, maxPos)]
    
    # Subset the array
    return(a[np.ix_(*n)])


def estimateQuantile(a, mp, nQuant):
    
    # Get all the seasons/years for the pixel of interest
    ref = a[mp[0], mp[1], :]
    ref = np.reshape(ref, (1, 1, len(ref)))
    
    am = array2TwoDim(a)
    
    i = 0
    
    while (np.sum(~np.isnan(ref)) < nQuant and i < max(a.shape[0:2])):
        i += 1
        ref = arrayAroundArray(a, mp, [i, i, a.shape[2]])
    
    if np.all(np.isnan(ref)):
        return(False)
    
    ref = array2TwoDim(ref)
    tmp = np.empty(ref.shape)
    tmp[:] = np.nan
    
    # Get the E c.d.f.
    for it in range(am.shape[1]):
        if not np.all(np.isnan(am[:, it])):
            tmp[:, it] = ECDF(am[:, it])(ref[:, it])
    
    tt = np.nanmean(tmp, axis=0)
    
    mu = np.nanmean(tt) 
    
    return(mu)


def gapFill(rasters, seasons, years, outFolder, suffix, nodata=None,
            iMax=np.inf,
            subsetSeasons=None, subsetYears=None, subsetMissing=None,
            clipRange=(-np.inf, np.inf)):
    '''
    Function to fill missing values in a series of rasters
    
    rasters (list): list of raster addresses
    seasons (list): Same size as rasters, gives for each the season. 
        Seasons should be numeric and in order by year. 
    years (list): Same size as rasters, gives for each the year.
        Years should be numeric and in order.
    outFolder (string): Address of output directory
    suffix (string): suffix to append to the name of the rasters after filling
    nodata (single num or list): Value to consider as no data in the rasters. 
        If None, will use the no data value in the raster if any
    fnSubset (function): function to use to subset the rasters
    fnPredict (function): function to use to predict the missing values. 
        The function should only return one element at a time
    iMax (int): Maximum number of iterations until missing is 
        returned as a predicted value
    subsetSeasons (list): list of seasons for which to predict missing values. 
        If None, all seasons are considered. Works as AND condition with 
        subsetYears
    subsetYears (list): list of years for which to predict missing values.
        If None, all years are considered. Works as AND condition with 
        subsetSeasons
    subsetMissing (logical array): Array of the same size as the rasters, 
        providing a subset of locations in which to consider missing values. 
        The same subset will be used for all seasons and all years.
    clipRange (tuple of length 2): Specifies the lower and the upper bound of the filled
        data. Values outside this range are clipped accordingly
        
    '''
    
    # Check the inputs
    
    # Check if all the raster addresses are correct and import them
    dst = []
    for r in rasters:
        if not os.path.isfile(r):
            return False
        
        try:
            dst.append(gdal.Open(r))
        except RuntimeError, e: 
            return False
        
        if dst[-1] is None:
            return False
    
    if not subsetYears: 
        subsetYears = set(years)
        subsetYears = list(subsetYears)
        subsetYears.sort()
    
    if not subsetSeasons:
        subsetSeasons = set(seasons)
        subsetSeasons = list(subsetSeasons)
        subsetSeasons.sort()
    
    # Import one band to get the array size
    a = dst[0].GetRasterBand(1).ReadAsArray()
    if not subsetMissing:
        # Import one band to get the array size
        subsetMissing = np.ones(a.shape, dtype=np.int)
    else:
        if not a.shape == subsetMissing.shape:
            return False
    
    if not nodata:
        nodata = dst[0].GetRasterBand(1).GetNoDataValue()
    
    if not nodata:
        return False
    
    if type(nodata) is int or nodata == np.nan:
        nodata = [nodata]
    
    p = 0
    # Loop through all the rasters in the subset to fill the missing values
    for y in subsetYears:
        for s in subsetSeasons:
            
            # Get the position of the raster for that year and season
            rIndex = zip(seasons, years).index((s, y))
            # Import that raster to get the positions of the missing values
            r = dst[rIndex].GetRasterBand(1).ReadAsArray()
            
            nodataOut = dst[rIndex].GetRasterBand(1).GetNoDataValue()
            if not nodataOut:
                nodataOut = nodata[0]
            
            # Get the position of all the missing data
            mps = []
            
            for n in nodata:
                miss = np.where(r == n)
                mps += map(list, zip(miss[0], miss[1]))
            
            # Loop through the missing data to replace in r
            for pixel in mps:
                
                # Create full position of pixel to replace
                mp = pixel + [s, y]
                
                i = 0 
                
                # Subset the rasters around the pixels
                if not nodataOut in nodata:
                    replaceVal = nodata + [nodataOut]
                else:
                    replaceVal = nodata
                    
                a = gapSubset(rasters=dst, seasons=seasons, years=years,
                             mp=mp, i=i, initialSize=[10, 10, 1, 5],
                             nodata=replaceVal)
                
                # Predict the value
                z = gapPredict(a=a[0], i=i, mp=a[1], nodataIn=nodata[0],
                               nodataOut=nodataOut)
                
                while z == nodataOut and i < iMax:
                    i += 1
                    
                    aNew = gapSubset(rasters=dst, seasons=seasons, years=years,
                                    mp=mp, i=i, initialSize=[10, 10, 1, 5],
                                    nodata=replaceVal)
                    if aNew[0].shape == a[0].shape:
                        break
                    
                    a = None
                    a = aNew
                    z = gapPredict(a=a[0], i=i, mp=a[1], nodataIn=nodata[0],
                               nodataOut=nodataOut)
                
                # Clip the value if needed
                if not z == nodataOut:
                    z = max(min(z, clipRange[1]), clipRange[0])
                
                # Replace the fitted value in the raster
                r[pixel[0], pixel[1]] = z
            
            # Prepare name for output
            outName = os.path.join(outFolder,
                                   re.sub('.tif', '_' + suffix + '.tif',
                                         os.path.basename(rasters[rIndex])))
            
            if os.path.isfile(outName):
                os.remove(outName)
             
            outR = md.new_raster_from_base(dst[rIndex], outName,
                                      'GTiff', nodataOut, gdal.GDT_Float32)
            
            outR.GetRasterBand(1).WriteArray(r)
            
            outR.FlushCache()
            outR = None


def gapPredict(a, i, mp, nodataIn, nodataOut=None, nTargetImage=5, nImages=4,
               nQuant=2):
    '''
    Function to predict a missing value
    
    a (array): Numeric array with three dimensions. Return value of gapSubset.
    i (int): The number of tried subsets that lead to a missing return
        value from gapPredict
    mp (list): list of 3 integers encoding the position of the missing 
        value in data to predict.
    nodataIn (num): value to be considered as no data in input data
    nodataOut (num): value to be used as no data for the output. 
        If None nodataIn will be used
    nTargetImage (int): Minimum number of non missing values in the image
        containing the missing value. If the criterion is not met,
        nodataOut is returned.
    nImages (int): Minimum number of non-empty images. If the criterion 
        is not met, nodataOut is returned.
    nQuant (int): Parameter passed to estimateQuantile
    '''
    
    outType = a.dtype
    
    a = a.astype(float)
    a[a == nodataIn] = np.nan
    
    am = array2TwoDim(a)
    
    if (np.sum(~np.isnan(am[:, mp[2]])) < nTargetImage):
        return(nodataOut)
    if (np.sum(np.apply_along_axis(np.any, 0, ~np.isnan(am))) < nImages): 
        return(nodataOut)
    
    tau = estimateQuantile(a=a, mp=mp, nQuant=nQuant)
    if not tau or tau == np.nan: 
      return(nodataOut)
    
    s = gapScore(mat=am)
    
    s = np.array(s)
    r = np.copy(s)
    
    r[~np.isnan(s)] = ss.rankdata(s[~np.isnan(s)])
    
    # Create a pandas dataframe to ru the quantile regression
    df = pd.DataFrame(data=np.column_stack([am.flatten('F'), np.repeat(r, repeats=am.shape[0])]),
                      columns=["z", "rank"])
    
    # Run quantile regression
    mod = smf.quantreg('z ~ rank', df) 
    
    m = mod.fit(q=tau, max_iter=800, p_tol=1e-02) 
    
    # Estimate fitted value
    p = m.params['Intercept'] + m.params['rank'] * r[mp[2]]
    
    if np.issubdtype(outType, np.integer):
        p = int(p)
    
    return(p)


def gapScore(mat):
    '''
    Compute the score for each column of the mat array
    mat array should have two dimensions only
    '''
    
    cols = mat.shape[1]
    rows = mat.shape[0]
    
    out = []
    for c in range(cols):
        # Extract column c
        colC = mat[:, c].reshape([rows, 1])
        colC = np.tile(colC, cols - 1)
        
        # Compare to the other columns
        o = np.divide(colC, np.delete(mat, c, 1))
        o[~np.isnan(o)] = o[~np.isnan(o)] > 1
        
        # Compute the mean
        out.append(np.nanmean(np.nanmean(o, 0)))
        
    return(out)


def gapSubset(rasters, seasons, years, mp, i, initialSize=[10, 10, 1, 5],
              nodata=None):
    out = arrayAroundRasters(rasters=rasters, seasons=seasons, years=years,
                mp=mp, size=[sum(x) for x in zip(initialSize, [i, i, 0, 0])],
                nodata=nodata)
    
    return out

    
if __name__ == '__main__':
    inDir = '/home/olivierp/jde_coffee/MODIS/collection6/Vietnam/LD/masked_missing'
    
    inputRasters = [os.path.join(inDir, f) for f in os.listdir(inDir) if
                    f.endswith('.tif')]
    inputRasters.sort()
    
    datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in inputRasters]
    # Transform into date format
    datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
    
    # Transform into days from start of the year
    days = [int(d.strftime('%j')) for d in datesAll]
    
    # Get the years for the files on disk
    years = [int(d.strftime('%Y')) for d in datesAll] 
    
    gapFill(rasters=inputRasters, seasons=days, years=years,
            outFolder='/home/olivierp/jde_coffee/MODIS/collection6/Vietnam/LD/filled_missing',
            suffix='f', nodata=[-3000], iMax=np.inf,
            subsetSeasons=None, subsetYears=None, subsetMissing=None,
            clipRange=(-2000, 10000))
    
    '''
    seasons = days
    
    dst = []
    for r in inputRasters:
        dst.append(gdal.Open(r))
        
    subsetYears = [2017]
    
    subsetSeasons = [337]
    
    # Import one band to get the array size
    a = dst[0].GetRasterBand(1).ReadAsArray()
    subsetMissing = np.ones(a.shape, dtype=np.int)
    a = None
    
    nodata = [-3000]
    
    y = subsetYears[0]
    s = subsetSeasons[0]
    
    # Get the position of the raster for that year and season
    rIndex = zip(seasons, years).index((s, y))
    # Import that raster to get the positions of the missing values
    r = dst[rIndex].GetRasterBand(1).ReadAsArray()
            
    nodataOut = dst[rIndex].GetRasterBand(1).GetNoDataValue()
    
    # Get the position of all the missing data
    mps = []
    for n in nodata:
        miss = np.where(r == n)
        mps += map(list, zip(miss[0], miss[1]))
    
    pixel = mps[3000]
    
    mp = pixel + [s, y]
    print(mp)  
    
    i = 0 
                
    # Subset the rasters around the pixels
    if not nodataOut in nodata:
        replaceVal = nodata + [nodataOut]
    else:
        replaceVal = nodata
                    
    A = gapSubset(rasters=dst, seasons=seasons, years=years,
                  mp=mp, i=i, initialSize=[10, 10, 1, 5],
                  nodata=replaceVal)
    
    z = gapPredict(a=A[0], i=i, mp=A[1], nodataIn=nodata[0],
                               nodataOut=nodataOut)
    
    # with open('/home/olivierp/jde_coffee/Temp/test.txt', 'wb') as f:
    #    np.savetxt(f, a[0], fmt='%.5f')
    '''
