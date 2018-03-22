# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Script to download MODIS images from the server, 
    and mosaic the images and mask them for each region
    to process.
    
    The image treatment process for each AOI is the following:
    1- Download the MODIS tiles
    2- Mosaic and cut the the extent of each AOI
    3- Mask missing/bad quality pixels in entire extent and mark those that 
        need to be filled
    4- Fill missing values
    5- Mask the image to only keep the pixels of interest
    6- Smooth the images temporally
    
"""

import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

import os, re, multiprocessing
import pathos.multiprocessing as mp
import MODIS_gedata_toolbox as md  # GEDATA toolbox for MODIS related tools
import gapfill_test  # Python implementation of the interpolation algorithm
from datetime import datetime
from csv import DictWriter
import time


def main():
    #####################################################################################################################
    #####################################################################################################################
    # #PARAMETERS
    
    # Allow parallel computing?
    allowPara = False
    # Number of cores to use?
    nCores = 3
    
    # Satellite
    satelliteModis = 'terra'  # 'terra' # 'aqua'
    
    # Root folder
    prefixRootSys = '/media/olivier/olivier_ext1/gedata_current/jde_coffee'  #'E:/gedata_current' #'/home/olivierp/jde_coffee'
    
    # #DIRECTORIES parameters
    # Working directory
    dst = os.path.join(prefixRootSys, 'MODIS/collection6/' + satelliteModis + '/Brazil')
    
    # #REGIONS parameters
    # Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    # Names of the regions (also folders names) 
    states = ["CER"]  # ["CER", "CHA", "CO", "ES", "MO", "SDM", "SP", "ZM"]
    # Name of the subfolder where to save the masked images
    statesMaskedFolder = 'masked_missing'
    # Name of the subfolder where to save the filled images
    statesFilledFolder = 'test_speed_filling'
    
    ############ FILL MISSING
    # Input folder for the images to fill
    inMissing = statesMaskedFolder
    # Output folder for the images to fill
    outMissing = statesFilledFolder
    # Year(s) of images to fill
    yearsMissing = [2018]
    # Day(s) of images to fill
    daysMissing = [[33, 49]]
    # Suffix to put at the end of the name of the 
    # images after filling
    suffMissing = 'f'
    # !!! The two conditions are additive (AND)    
    
    #####################################################################################################################
    #####################################################################################################################
    # #ACTIVE CODE
    
    if allowPara and not nCores:
        nCores = multiprocessing.cpu_count()
    
    print('Starting interpolation of missing values')
    
    for s in states:
        print('   Starting region ' + s)
        
        inputRasters = [os.path.join(dst, s, inMissing, f) for 
                        f in os.listdir(os.path.join(dst, s, inMissing)) 
                        if f.endswith('.tif')]
        
        inputRasters.sort()
        
        # Get the dates
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in inputRasters]
        # Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
        
        # Transform into days from start of the year
        days = [int(d.strftime('%j')) for d in datesAll]
        
        # Get the years for the files on disk
        years = [int(d.strftime('%Y')) for d in datesAll] 
        t0 = time.time()
        expans = gapfill_test.gapFill(rasters=inputRasters, seasons=days, years=years,
                                 outFolder=os.path.join(dst, s, outMissing),
                                 suffix=suffMissing, nodata=[-3000], iMax=20,
                                 subsetSeasons=daysMissing, subsetYears=yearsMissing,
                                 subsetMissing=None, clipRange=(-2000, 10000),
                                 parallel=allowPara, nCores=nCores)
        print time.time() - t0
        # print expans


if __name__ == '__main__':
    main()
