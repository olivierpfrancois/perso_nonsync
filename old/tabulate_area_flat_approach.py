# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

"""
The algorithm loops through each pixel of the raster, compute the coordinates of
the centroid of the pixel, looks which polygons the pixel falls into,
and adds its value to each of the polygons it intersects.
The tabulated values are exported into a tab delimited file.

The resulting algorithm is very slow.

"""

# Import the libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
import sys, os, processing, numpy as np
from osgeo import ogr, gdal, osr
from collections import Counter #Tool to get unique value counts from list
from csv import DictWriter


def run_script(iface):
	
	
	### GET PARAMETERS
	
    #Define the folder address with the shapefile
	shpFolder = 'D:\TabulateAreaTEst' #'/media/olivier/olivier_ext/TabulateAreaTEst' # args['pathshp']
	#Define the folder address with the image
	imgFolder = 'D:\TabulateAreaTEst' #'/media/olivier/olivier_ext/TabulateAreaTEst' #args['pathtif']
	
	#Name of the shapefile
	shpName = 'grid_IA.shp' #args['shp']
	#Name of the image
	imgName = 'CDL2015_IA.tif' #args['tif']
	
	##No data value
	#if 'nodata' in args:
	#	no_data = args['nodata']
	#else:
	#	no_data = set()
	no_data = set([0])
	
	##Index to identify the polygons
	#if 'index' in args:
	#	index = args['index']
	#else:
	#	index = None
	index = 'GridCellID'
	
	##Prefix for the data
	prefix = 'cat'
	
	############################
	
	### RUN ANALYSIS
	
	#Import the raster
	ds = gdal.Open(imgFolder+'/'+imgName)
	
	#Get the no data value and add with user-supplied if not the same
	if ds.GetRasterBand(1).GetNoDataValue():
		no_data.add(ds.GetRasterBand(1).GetNoDataValue())
	no_data = list(no_data)
	
	#Import the shapefile
	driver = ogr.GetDriverByName('ESRI Shapefile')
	dataSource = driver.Open(shpFolder+'/'+shpName, 0) # 0 means read-only. 1 means writeable.
	# Check to see if shapefile is found.
	if dataSource is None:
		print 'Could not open %s' % (shpFolder+'/'+shpName)
	else:
		layer = dataSource.GetLayer(0)
	
	if not index == None:
		#Get the index of the field with the polygon ids
		idx = layer.GetLayerDefn().GetFieldIndex(index)
	else:
		#Create an index within the polygon layer
		#Get the index for that new variable
		pass
	
	#Crop the raster by the layer extent in case it is smaller
	args = {"INPUT":imgFolder+'/'+imgName, "PROJWIN":", ".join(map(str, layer.GetExtent())), 
					"RTYPE":4, "COMPRESS":0, "BIGTIFF":3, "OUTPUT":None}
	if len(no_data) > 0:
		args["NO_DATA"] = no_data[0]
	
	ds = processing.runalg("gdalogr:cliprasterbyextent", args)
	ds = ds['OUTPUT'] #Get the address of the temporary file with the result of the algorithm
	ds = gdal.Open(ds) #Import again
	
	#Extract the raster georeference information
	geo_transform = ds.GetGeoTransform()
	xOrigin = geo_transform[0] #Top left x
	yOrigin = geo_transform[3] #Top left y
	xRotation = geo_transform[2] #Rotation, 0 if image is north up
	yRotation = geo_transform[4] #Rotation, 0 if image is north up
	pixelWidth = geo_transform[1] #W-E pixel resolution
	pixelHeight = geo_transform[5] #N-S pixel resolution

	#Get raster projection
	proj_ds = ds.GetProjection()	
	
	#Get the data type for the array
	d_type = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
	
	#Get image size
	rows = ds.RasterYSize
	cols = ds.RasterXSize
	
	#Get the first band
	ds_band = ds.GetRasterBand(1)
	
	#Prepare an empty dictionary to hold the results
	tabulate = {}
	
	#Loop through the rows of the raster and process all pixels individually
	for L in range(rows):
		data = ds_band.ReadAsArray(0, L, cols, 1)
		for P in range(data.shape[1]):
			#Get the approximate coordinates of the center of the pixel
			Xp = xOrigin + P*pixelWidth + L*xRotation + pixelWidth/2
			Yp = yOrigin + P*yRotation + L*pixelHeight + pixelHeight/2
		
			# DEBUG: take a look at the coords..
			#print "Coords: " + str(x) + ", " + str(y)
			
			#Transform into a point geometry
			pt = ogr.Geometry(ogr.wkbPoint)
			pt.SetPoint_2D(0, Xp, Yp)
			
			#Set a filter to only see the polygons that intersect
			layer.SetSpatialFilter(pt)
			
			#Loop through intersected polygons to add the value of the raster pixel
			for feat_in in layer:
				#Get the id
				id_in = feat_in.GetFieldAsString(idx)
				
				try:
					#Add to the result dictionary if the key already exists
					#if str(data[L]) in tabulate[id_in]:
					tabulate[id_in][str(data[0,P])] += 1
				except:
					#Create the key if does not exists
					#else:
					try:
						tabulate[id_in][str(data[0,P])] = 1
					except:
						#Add the key to the dictionary and prepare it as an empty dictionary
						tabulate[id_in] = {}
						tabulate[id_in][str(data[0,P])] = 1
	
	
	#Transform the result dictionary for the export
	out = []
	#Prepare a set to hold the unique values in the raster
	uniques = set()
	for k, v in tabulate.items():
		v['id'] = k
		uniques.update(v.keys())
		out.append(v)
	
	#Export to a tab delimited file
	with open(shpFolder+'/test_flat.txt', "w") as f:
		dict_writer = DictWriter(f, ['id']+sorted(uniques), delimiter="\t")
		dict_writer.writeheader()
		for p in out:
			dict_writer.writerow(p)
	
	'''
	#Remove any layers with same name in registry if already one
	for layer in QgsMapLayerRegistry.instance().mapLayers().values():
		if layer.name() == "Polygons" or layer.name() == "Image":
			QgsMapLayerRegistry.instance().removeMapLayers([layer.id()])
					
	#Import the ground shapefile of polygons
	grid = iface.addVectorLayer(shpFolder+"/"+shpName, "Polygons", "ogr")
	if not grid or not grid.isValid():
			print("Failed to load %s" %(shpFolder+"/"+shpName))
			sys.exit()
            
	#Import the raster image
	img = iface.addRasterLayer(imgFolder+"/"+imgName, "Image")
	if not img or not img.isValid():
			print("Failed to load %s" %(imgFolder+"/"+imgName))
			sys.exit()
	
	#Add an index to be able to identify the vt polygons uniquely
	#if no index was specified by the user.
	i = 1
	
	# Create a spatial index for the shapefile features to increase the computation speed
	# Get all the features of the cell layer
	gridFeatures = {feature.id(): feature for (feature) in grid_shp.getFeatures()}
	# Build the spatial index for faster lookup.
	gIndex = QgsSpatialIndex()
	map(gIndex.insertFeature, gridFeatures.values())
	
	#Loop over all the raster cells to extract the cell id it belongs to and update a count
	#get the raster data
	prov = img.dataProvider()
	width = prov.xSize()
	height = prov.ySize()
	block = prov.block(1, prov.extent(), width, height)
	#Add the no data value from the raster if any
	if block.noDataValue():
		no_data.add(block.noDataValue())
	#do the loop over all pixels
	for x in range(width):
		for y in range(height):
			
			Xp = padfTransform[0] + P*padfTransform[1] + L*padfTransform[2];
			Yp = padfTransform[3] + P*padfTransform[4] + L*padfTransform[5];

			values.append(int(block.value(x, y)))
	'''
