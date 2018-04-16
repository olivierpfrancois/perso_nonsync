# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

"""
The algorithm imports the raster and transforms it into a numpy 
array. It loops through the polygons in the vector layer, crops the 
raster to the polygon, then loops through all the pixels of the cropped
raster to count those that fall into the polygon.

The tabulated values are exported into a tab delimited file.
"""

# Import the libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
import sys, os, processing, numpy as np
from osgeo import ogr, gdal, osr
import collections #Tool to get unique value counts from list
from csv import DictWriter
import shapely
from shapely import wkt
import time
def run_script(iface):
	
	def GetExtentRaster(gt,cols,rows):
		''' Return list of corner coordinates from a geotransform
			@type gt:   C{tuple/list}
			@param gt: geotransform
			@type cols:   C{int}
			@param cols: number of columns in the dataset
			@type rows:   C{int}
			@param rows: number of rows in the dataset
			@rtype:    C{[float,...,float]}
			@return:   coordinates of each corner
		'''
		ext=[]
		xarr=[0,cols]
		yarr=[0,rows]
		
		#Goes through each of the four corners to compute the coordinates
		for px in xarr:
			for py in yarr:
				x=gt[0]+(px*gt[1])+(py*gt[2])
				y=gt[3]+(px*gt[4])+(py*gt[5])
				ext.append([x,y])
			yarr.reverse()
		return ext
	
	def world2Pixel(gt, x, y):
		"""
		Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
		the pixel location of a geospatial coordinate
		
		The returned pixel may be negative if the coordinate is 
		outside of the extent of the geotransform to the north or
		west.
		
		Starts from the upper left corner, so the line value increases going down.
		"""
		ulX = gt[0]
		ulY = gt[3]
		xDist = gt[1]
		yDist = gt[5]
		pixel = int((x - ulX) / xDist)
		line = int((ulY - y) / xDist)
		return (pixel, line)
	
	
	############################
	##### Active code
	
	
	### GET PARAMETERS
	t0 = time.time()
	#Define the folder address with the shapefile
	shpFolder = '/media/olivier/olivier_ext/TabulateAreaTEst' # args['pathshp'] #'D:\TabulateAreaTEst'
	#Define the folder address with the image
	imgFolder =  '/media/olivier/olivier_ext/TabulateAreaTEst' #args['pathtif'] #'D:\TabulateAreaTEst'
	
	#Name of the shapefile
	shpName = 'grid_IA.shp' #'smallpols.shp' #args['shp']
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
	index = 'GridCellID' #'ID'
	
	##Prefix for the data
	prefix = 'cat_'
	
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
	
	#Get the extent of the shapefile layer
	shp_ext = layer.GetExtent()
	
	# Get Raster extent coordinates
	rast_ext = GetExtentRaster(ds.GetGeoTransform(), ds.RasterYSize, ds.RasterXSize)
	# Format to the same tuple structure as the vector layer
	rast_ext = (rast_ext[0][0],rast_ext[2][0],rast_ext[1][1],rast_ext[0][1])
	
	#Compare the two extent to see if raster inside the shapefile
	if (rast_ext[0] < shp_ext[0] or rast_ext[1] > shp_ext[1] or 
		rast_ext[2] < shp_ext[2] or rast_ext[3] > shp_ext[3]):
		
		#Crop the raster by the layer extent in case it is smaller
		args = {"INPUT":imgFolder+'/'+imgName, "PROJWIN":", ".join(map(str, layer.GetExtent())), 
					"RTYPE":4, "COMPRESS":0, "BIGTIFF":3, "OUTPUT":None}
		if len(no_data) > 0:
			args["NO_DATA"] = no_data[0]
		
		ds = processing.runalg("gdalogr:cliprasterbyextent", args)
		ds = ds['OUTPUT'] #Get the address of the temporary file with the result of the algorithm
		ds = gdal.Open(ds) #Import again
	
	#Extract the raster projection information
	geo_transform = ds.GetGeoTransform()
	proj_ds = ds.GetProjection()
	xOrigin = geo_transform[0] #Top left x
	yOrigin = geo_transform[3] #Top left y
	xRotation = geo_transform[2] #Rotation, 0 if image is north up
	yRotation = geo_transform[4] #Rotation, 0 if image is north up
	pixelWidth = geo_transform[1] #W-E pixel resolution
	pixelHeight = geo_transform[5] #N-S pixel resolution
	
	d_type = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
	
	#Convert the first band to an array
	try:
		ds_band = np.array(ds.GetRasterBand(1).ReadAsArray())
		if len(ds_band.shape) == 3:
			ds_band = ds_band[:,:,0] #Keep only the first two dimensions
	except ValueError:
		print("Raster file is to big for processing. Please crop the file and try again.")
		return
	
	#Close the raster
	ds = None
	
	##Crop the raster to the extent of the layer in case it is smaller
	#ds_band = cropRasterArrayToLayer(ds_band, layer.GetExtent(), geo_transform)
	#if not ds_band:
	#	print("The layer extent and the raster don't overlap.")
	#	return
	
	#Get number of rows, columns and bands of the raster
	rows, cols = ds_band.shape
	
	# Get Raster extent coordinates
	extent = GetExtentRaster(geo_transform,rows,cols)
	# Format to the same tuple structure as the vector layer
	extent = (extent[0][0],extent[2][0],extent[1][1],extent[0][1])
	
	#Prepare a dictionary to hold the counts
	tabulate = {}
	
	#Iterate over the features
	for feat in layer:
		
		# Get geometry of feature
		geom = feat.GetGeometryRef()
		
		#Get the id
		id_in = feat.GetFieldAsString(idx)
		
		#Get geometry envelope
		geom_extent = geom.GetEnvelope()
		
		#Transform into pixel and line positions
		pixels = []
		for ptsX in geom_extent[0:2]:
			for ptsY in geom_extent[2:]:
				pixels.append(world2Pixel(geo_transform, ptsX, ptsY))
		
		#Get pixel extent
		xpix, ypix = zip(*pixels)
		pixExtent = [min(xpix), max(xpix)+1, min(ypix), max(ypix)+1]
		
		#Remove negative values
		pixExtent = [max(0,i) for i in pixExtent]
		pixExtent[0] = min(pixExtent[0],cols)
		pixExtent[1] = min(pixExtent[1],cols)
		pixExtent[2] = min(pixExtent[2],rows)
		pixExtent[3] = min(pixExtent[3],rows)
		
		if ((pixExtent[0] == rows+1 and pixExtent[2] == cols+1) or 
			(pixExtent[1] == 0 and pixExtent[3] == 0)):
			#Polygon is outside of raster
			continue
		
		#Transform polygon into shapely object
		poly = wkt.loads(geom.ExportToWkt())
		
		#Loop through all the raster pixels in that extent, 
		#Get the coordinates and the values
		pts = []
		val = []
		for L in range(pixExtent[2],pixExtent[3]):
			for P in range(pixExtent[0],pixExtent[1]):
				#Get coordinates of the center of the pixel
				Xp = xOrigin + P*pixelWidth + L*xRotation + pixelWidth/2
				Yp = yOrigin + P*yRotation + L*pixelHeight + pixelHeight/2
				
				pts.append((Xp,Yp))
				val.append(str(ds_band[L,P]))
		
		#Transform into shapely points
		pts = [shapely.geometry.Point(i,j) for i,j in pts]
		#print("nb points: "+str(len(pts)))
		#Get the indexes of the points that fall into the polygon
		inpts =  [i for i, pt in enumerate(pts) if poly.intersects(pt)] #pt.within(poly)
		#print("index: "+' '.join(map(str,inpts)))
		#Get the values for those points and add prefix
		values = [prefix+val[i] for i in inpts]
		
		#Extract the unique values and counts
		tabulate[id_in] = collections.Counter(values)
	
	#Transform the result dictionary for the export
	out = []
	#Prepare a set to hold the unique values in the raster
	uniques = set()
	for k, v in tabulate.items():
		uniques.update(v.keys())
		v['id'] = k
		out.append(v)
	
	#Export to a tab delimited file
	with open(shpFolder+'/test_flat.txt', "w") as f:
		dict_writer = DictWriter(f, ['id']+sorted(uniques, key = lambda x: int(x.split("_")[1])), delimiter="\t", restval="0")
		dict_writer.writeheader()
		for p in out:
			dict_writer.writerow(p)
	t1 = time.time()
	print(str(t1-t0))
	'''
	To do:
	-Cut in pieces automatically if too big to process, or work by block...
	-Check if raster has projection
	-Check if shapefile has projection
	-Check if projections are the same
	'''
	
	
	#Get the epsg code of the layer projection
	#crs = layer.GetSpatialRef()
    #crs.AutoIdentifyEPSG()
    #code = srs.GetAuthorityCode(None)
    #if code:
    #    srid= code
	
	
	'''
	#If the index is none in the shapefile, create one
	if index == None:
		#TODO
		pass
		
	#Rasterize the shapefile
	driver = gdal.GetDriverByName('MEM')  # In memory dataset
    mask = driver.Create('', cols, rows, 1, gdal.GDT_UInt16)
    mask.SetGeoTransform(geo_transform)
    mask.SetProjection(projection)
    gdal.RasterizeLayer(mask, [1], layer, options = ["ATTRIBUTE="+index])
    
	mask_band = np.array(mask.GetRasterBand(1).ReadAsArray())
	mask = None #Close the raster to free memory space
	
	#Make sure the rasterized shapefile has the same extent, columns and rows as the image
	
	#Loop through the array indices[] to export the values and polygon indices
	
	#Export the data
	'''
