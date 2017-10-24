# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" 
Script to rasterize a regular grid of square cells.
It only works for utm projections.
The three inputs are the 
-full address of a regular grid [text]
-name of the variable to use for the raster values [text]
-resolution of the grid in meters [number]
-resolution of the raster in meters [number]

The code takes the grid, transforms it into a point shapefile
using the coordinates of the centroid for each cell, 
rasterizes the point shapefile at the resolution of the grid,
and changes the resolution of the resulting file to the wanted
resolution for the raster.

This is a workaround for gdal_rasterize, which otherwise 
creates a raster that does not perfectly overlaps with the grid.

 """

# Imports

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
import processing, re, tempfile, os

def run_script(iface):
    
	#INPUTS:
	nm = "/home/olivier/Desktop/temp/AOI_grid_ZM_2km.shp"
	nmVar = "serp_id"
	resGrid = 2000
	resRaster = 6
	
	#FUNCTIONS
	def toRaster(nm, nmVar, resGrid, resRaster):
		
		#Import the grid
		grid = QgsVectorLayer(nm, "region", "ogr")
		
		#Get the fields names
		field_names = [field.name() for field in grid.pendingFields()]
		
		if not nmVar in field_names:
			print("Variable name provided is not one of the grid's fields.")
			return
		
		#Create the name of the output raster
		outName = re.sub(".shp",".tif",nm)
		
		#The gdal rasterize function works from the corners of the cells
		#and the resulting raster is off center by half a pixel.
		#So I have to first create a points shapefile using the 
		#centroids, then rasterize at the same resolution as the grid,
		#and finally change the resolution to the final resolution.
		
		#Create point layer from grid centroids
		pointsLayer = QgsVectorLayer("Point?crs=EPSG:"+str(grid.crs().postgisSrid()), "gridPoints", "memory")
		pr = pointsLayer.dataProvider()
		#Enter editing mode
		pointsLayer.startEditing()
		#Add serp_id field
		pr.addAttributes([QgsField(nmVar, QVariant.Int)])
		#Add features
		#Loop over all cells in the grid
		fets = []
		i = 0
		for cell in grid.getFeatures():
			#Create an empty feature
			fet = QgsFeature()
			#Declare it as point and add the coordinates
			fet.setGeometry(cell.geometry().centroid())
			#Prepare the attributes and add them
			fetAttr = [cell[nmVar]]
			fet.setAttributes(fetAttr)
			#Add to list
			fets.append(fet)
		#Add the features
		pr.addFeatures(fets)
		#Commit changes
		pointsLayer.commitChanges()
		
		extent = pointsLayer.extent()
		extent = (extent.xMinimum(), extent.xMaximum(), extent.yMinimum(), extent.yMaximum())
		
		#Save to a temporary file
		tempout = str(os.path.join(tempfile.gettempdir(),"temp.shp"))
		QgsVectorFileWriter.writeAsVectorFormat(pointsLayer, tempout, "utf-8", None, "ESRI Shapefile")
		#QgsMapLayerRegistry.instance().addMapLayer(pointsLayer)
		
		first = processing.runalg("gdalogr:rasterize", {"INPUT":tempout, "FIELD":nmVar, 
				"DIMENSIONS":1, "WIDTH":resGrid, "HEIGHT":resGrid, "TFW":False, "RTYPE":2, 
				"NO_DATA":0, "BIGTIFF": 2, "OUTPUT":None})
		first = first['OUTPUT']
		
		#Remove temporary file
		#if os.path.exists(tempout):
		#	ogr.GetDriverByName("ESRI Shapefile").DeleteDataSource(tempout)
		
		#Remove output file if it already exists
		try:
			os.remove(outName)
		except OSError:
			pass
		
		processing.runalg("gdalogr:warpreproject", {"INPUT":first, 
				"DEST_SRS":"EPSG:"+str(grid.crs().postgisSrid()), "NO_DATA":0, 
				"TR":resRaster, "METHOD":0, "RTYPE":2, "BIGTIFF": 2, 
				"OUTPUT":outName})
	
	#WORKING BIT
	toRaster(nm, nmVar, resGrid, resRaster)
