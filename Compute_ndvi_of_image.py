# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
This script uses the QGIS Raster Calculator to calculate the ndvi

This script takes four arguments, 
-the folder address of a tif image,			path=[string]
-the name of a tif image (with extension), 	tif=[string]
-the number for the IR band and, 			ir=[int]
-the number for the red band. 				red=[int]

It returns an image with a single band that is the ndvi.

The name of the resulting image is the same as the image, with a 
"_ndvi" suffix.
"""


# Some commonly used imports

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from qgis.analysis import *
import sys

def run_script(iface, **args):
	
	#Define the folder address with the image
	imgFolder = args['path']
	
	#Name of the image
	imgName = args['tif']
	
	#Bands for the calculation
	irB = args['ir']
	redB = args['red']
	
	#Remove any layers with same name in registry if already one
	for layer in QgsMapLayerRegistry.instance().mapLayers().values():
		if layer.name() == "Image":
			QgsMapLayerRegistry.instance().removeMapLayers([layer.id()])
	
	#Import the raster image
	img = iface.addRasterLayer(imgFolder+"/"+imgName, "Image")
	imgName = imgName.replace(".TIF", ".tif")
	
	if not img or not img.isValid():
			print("Failed to load %s" %(args['path']+'/'+args['tif']))
			sys.exit()
	
	###Compute the NDVI
	print("Starting the computation of the NDVI")
	
	###Using the gdal raster calculator
	#img = processing.runalg("gdalogr:rastercalculator", 
	#			{"INPUT_A":img, "BAND_A":str(irB), 
	#			"INPUT_B":img, "BAND_B":str(redB),  
	#			"FORMULA":"(A-B)/(A+B)",  
	#			"RTYPE":5, "OUTPUT":None})
	#img = img['OUTPUT']
	#img = QgsRasterLayer(img, "Image")
	
	##Using the QGIS raster calculator
	# create raster calculator objects for the red and IR bands
	ir = QgsRasterCalculatorEntry()
	red = QgsRasterCalculatorEntry()
	# connect the raster calculator objects to the raster (same here for both)
	ir.raster  = img
	red.raster = img
	# tell the raster calculator bands to use
	ir.bandNumber  = irB
	red.bandNumber = redB
    # construct input image variable names for raster calculator ("name@band")
	ir.ref  = "Image" + '@' + str(irB)
	red.ref = "Image" + '@' + str(redB)
    # construct raster calculator formula
    # (multiply numerator by 1.0 to force it to be a real number)
	ndvi_formula = ('(1.0 * (%s - %s)) / (1.0 * (%s + %s))' % 
				(ir.ref, red.ref, ir.ref, red.ref))
    # assemble a calculation description
	ndvi_calc = QgsRasterCalculator(ndvi_formula, imgFolder+imgName.replace(".tif","_ndvi.tif"), "GTiff",
                               img.extent(), img.width(), img.height(), [ir, red])
    # run the NDVI calculation
	ndvi_calc.processCalculation()
	# import the result
	#img = QgsRasterLayer(imgFolder+imgName.replace(".tif","_ndvi.tif") "Image")
	
	#Remove layer from qgis
	for layer in QgsMapLayerRegistry.instance().mapLayers().values():
		if layer.name() == "Image":
			QgsMapLayerRegistry.instance().removeMapLayers([layer.id()])
	
	print('NDVI raster successfully created in %s' %(imgFolder+imgName.replace(".tif","_ndvi.tif")))
