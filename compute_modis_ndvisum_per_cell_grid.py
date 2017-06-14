# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
    Takes a batch of modis ndvi images in a folder
    and a shapefile grid, 
    and computes the sum of the ndvi for each grid cell.
    The output is exported in a tab delimited file.
    
    For a new grid it is alsoo necessary to count the number
    of modis pixels falling into each grid cell, so to 
    run the same script on a modis image where all the
    valid pixels have been set to 1.
    (In that case need to change code at line 78)
    
    The script is used for the analysis of the modis data 
    together with weather and production.
"""

# Some commonly used imports

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from osgeo import ogr, gdal
from tabulate_area_functions import *
import re, os
from csv import DictWriter

def run_script(iface):
    #shpToRaster('D:/gedata_current/jde_coffee/MODIS/ES/ES_areas_segments_modis_sum.shp', 'hectares', 2000)
    
    '''
    Need to also compute the MODIS weights for the grid (Nb of MODIS pixels falling into every grid)
    This is simply done by tabulating with the sum statistics a MODIS image where all the pixels have been set to 1.
    '''
    
    #################################################################################################################################
    #################################################################################################################################
    ### GET PARAMETERS
    #Address of the shapefile
    shpName = 'D:/gedata_current/jde_coffee/data/ES/areas/ES_direct_cells_conf_cells1K.shp'
    
    #Name of variable in the shapefile with the ID of the polygons
    index = 'serp_id' #'Station'
    
    #Define the folder address with the MODIS images from which to extract the data
    imgFolder = 'D:/gedata_current/jde_coffee/MODIS/collection6/ES/smooth_data'
    #imgFolder = 'D:/gedata_current/jde_coffee/MODIS/collection6/ES/weights' #Folder for the computation of MODIS weights
    
    #Name for the txt output with the sum of NDVI inside each polygon
    outName = 'D:/gedata_current/jde_coffee/MODIS/collection6/ES/weights/ES_1kmgrid_modis_weights_sum.txt'
    #outName = 'D:/gedata_current/jde_coffee/MODIS/collection6/ES/weights/ES_1kmgrid_modis_weights_sum.txt' #Name for the computation of MODIS weights
    
    
    
    #################################################################################################################################
    #################################################################################################################################
    ### ACTIVE CODE
    
    #Get all the files in the folder
    imgNames = os.listdir(imgFolder)
    imgNames = [img for img in imgNames if img.endswith('.tif')]
    
    tot = str(len(imgNames))
    
    #Names for the suffixes for the variables created by the tabulate area function
    
    imgSuffixes = []
    for nm in imgNames:
        suff = re.search('_[0-9]{2}([0-9]{2}-[0-9]{2}-[0-9]{2})_', nm).group(1)
        suff = suff.replace('-','')
        imgSuffixes.append(suff)
    
    #imgSuffixes = ['weight'] #Suffix for the computation of modis weights
    
    ##Option for the rasterization
    alltouch = 'FALSE'
    
    ### RUN ANALYSIS
    data = []
    
    i = 1
    for img in imgNames:
        print('Processing image '+str(i)+' of '+tot)
        i+=1
        
        #Extract the statistics for the grid cells
        out = tabulate_area(rasterName=imgFolder+'/'+img, bands=[1], shpName=shpName, 
                    index=index, tabulate=False, prefix='', statistics=['sum'],
                    outTxt=None, addShp=False, addSuffix='', alltouch=alltouch)
        #The output is a dictionary of which we only want the statistics part.
        #That part is a list with one element per band. Here there is only one.
        data.append(out['statistics'][0])
    
    #Combine the data for each of the outputs
    #data is a list of dictionaries with the same keys (one per polygon)
    d = {}
    for k in data[0].iterkeys():
        d[k] = [cell[k] for cell in data]
    #Each element of d is now a list of dictionaries with each one statistics
    out=[]
    for K,V in d.items():
        out.append({}) 
        #Combine all the statistics into one dictionary
        for dic,s in zip(V,imgSuffixes):
            out[-1].update({k+s:v for k,v in dic.items()})
        out[-1]['id'] = K
    
    uniques = out[0].keys()
    uniques.remove('id')
    uniques.sort()
    
    with open(outName, "w") as f:
        dict_writer = DictWriter(f, ['id']+uniques, extrasaction='ignore', delimiter="\t", restval="0")
        dict_writer.writeheader()
        for p in out:
            dict_writer.writerow(p)
    