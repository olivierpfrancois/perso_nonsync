# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" Your Description of the script goes here """

# Some commonly used imports

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from tabulate_area_functions import *
import re, os

def run_script(iface):
	
	### GET PARAMETERS
	#Define the folder address with the shapefile
	shpFolder = 'D:/TabulateAreaTEst/test_MO_stats_ndvi' # args['pathshp'] #'/media/olivier/olivier_ext/TabulateAreaTEst'
	#Define the folder address with the image
	imgFolder =  'D:/TabulateAreaTEst/test_MO_stats_ndvi' #args['pathtif'] #'/media/olivier/olivier_ext/TabulateAreaTEst'
	
	#Name of the shapefile
	shpName = 'VT_mo01_2.shp' #'smallpols.shp' #args['shp']
	#Name of the images
	imgNames = ['01MO_20160522.TIF','L220_073-074-075_2015-10-16.tif','L220_074-075_2016-07-14.tif'] #args['tif']
	#Names for the suffixes for the variables in the shapefile
	imgSuffixes = ['_0522','_1016','_0714']
	
	##Index to identify the polygons
	#if 'index' in args:
	#	index = args['index']
	#else:
	#	index = None
	index = 'vt_id' #'ID' 
	
	##Prefix for the data
	prefix = 'cat_'
	
	##Option for the rasterization
	alltouch = 'FALSE'
	
	### RUN ANALYSIS
	for img, s in zip(imgNames, imgSuffixes):
		#Compute the ndvi
		ndvi_calc(imgFolder+'/'+img, redBand=1, irBand=4, 
					outTif=True, outName=None)
		
		ndvi = imgFolder+'/'+re.sub("(.tif|.TIF)", "_ndvi.tif", img)
		
		#Extract the statistics for the VT polygons
		tabulate_area(rasterName=ndvi, bands=[1], shpName=shpFolder+'/'+shpName, 
						index=index, tabulate=False, prefix=prefix, statistics=['mean','sd'],
						outTxt=re.sub(".tif", ".txt", ndvi), addShp=True, addSuffix=s, alltouch=alltouch)
		
		os.remove(ndvi)
		#os.remove(ndvi+'.aux.xml')
		
		