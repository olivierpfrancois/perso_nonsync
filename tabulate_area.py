# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

"""
The script takes as arguments
-the folder address of a shp of polygons, 		pathshp=[string]
-the name of a shp of polygons (w/ extension),	shp=[string]
-the folder address of a tif,					pathtif=[string]
-the name of a tif (w/ extension), and 			tif=[string]
-the nodata value for the image (optional)		nodata=[number]
-the name of the variable holding the unique index to
	identify the polygons (optional)			index=[string]

It loops through each polygon and layer in the image, 
and returns a .csv with the tabulate for each polygon.
If no index was provided, it creates an index in the 
order of the polygons in the shapefile.
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
		"""
		ulX = gt[0]
		ulY = gt[3]
		xDist = gt[1]
		yDist = gt[5]
		pixel = int((x - ulX) / xDist)
		line = int((ulY - y) / xDist)
		return (pixel, line)
	
	def cropRasterArrayToExtent(srcArray, lyrExtent, gt):
		'''
		Returns the source raster (as a numpy array) cropped to
		the extent of the layer, if the layer extent is smaller, 
		using the geotransform of the raster.
		
		The layer extent should be a list as provided by .GetExtent(),
		do in the order minX, maxX, minY, maxY
		'''
		
		#Get the pixel locations of the extent of the layer
		minX, maxX, minY, maxY = lyrExtent
		ulX, ulY = world2Pixel(gt, minX, maxY)
		lrX, lrY = world2Pixel(gt, maxX, minY)
		
		#Replace pixel position by 0 if the position is negative
		ulX = max(ulX, 0)
		ulY = max(ulY, 0)
		lrX = max(lrX, 0)
		lrY = max(lrY, 0)
		
		#Get number of rows, columns and bands
		rows, cols, n_bands = srcArray.shape
		
		#If all positions are out of the array, return None
		if max(lrX,lrY) == 0 or (ulY >= rows and ulX >= cols):
			return 

		# Clip the raster to the shapes bounding box
		clip = srcArray[ulY:lrY, ulX:lrX]
		
		return (clip)
	
	def getPtGeom(geom):
		#Lists to recover the coordinates of the points
		pointsX = []; pointsY = []
		
		#Get the geometry for this polygon
		ring = geom.GetGeometryRef(0)
		#Get number of points
		numpoints = ring.GetPointCount()
		#Extract coordinates
		for p in range(numpoints):
			lon, lat, z = ring.GetPoint(p)
			pointsX.append(lon)
			pointsY.append(lat)	
	
	def imageToArray(i):
		"""
		Converts a Python Imaging Library array to a
		gdalnumeric image.
		"""
		try:
			a = np.fromstring(i.tostring(),'b')
		except AttributeError:
			try:
				a = np.fromstring(i.tobytes(),'b')   
			except SystemError:
				return
		if a != None:
			a.shape=i.im.size[1], i.im.size[0]
		return a
	
	def getClipArray(layer, poly, srcArray, gt):
		'''
		Returns an array from a given polygon feature
		Assumes that the polygon is inside the raster extent
		@param layer: layer to which the polygon belongs
		@param poly: polygon feature from the layer
		@srcArray: raster layer (as array from gdal) that will be clipped
		@param gt: gdal geomatrix (gdal.GetGeoTransform()) from underlying raster
		'''
		
		# Get the geometry ref from the polygon
		try:
			geom = poly.GetGeometryRef()
		except AttributeError:
			return None
		
		#Check if it is a multipolygon
		if geom.GetGeometryCount() > 1:
            # TODO: What to do with multipolygons?
            print("The Shapefile contains Multipolygons which are not supported")
			return None

		#Get the pixel locations of the extent of the layer
		minX, maxX, minY, maxY = layer.GetExtent()
		ulX, ulY = world2Pixel(gt, minX, maxY)
		lrX, lrY = world2Pixel(gt, maxX, minY)
		# Calculate the pixel size of the new image
		pxWidth = int(lrX - ulX)
		pxHeight = int(lrY - ulY)
		
		# If the clipping features extend out-of-bounds and ABOVE the raster...
		if gt[3] < maxY:
			# In such a case... ulY ends up being negative--can't have that!
			iY = ulY
			ulY = 0

		# Clip the raster to the shapes boundingbox
		clip = srcArray[ulY:lrY, ulX:lrX]
		
		# Create a new geomatrix for the image that is covered by the layer
		geoTrans = list(gt)
		geoTrans[0] = minX
		geoTrans[3] = maxY
		
		# Map points to pixels for drawing the boundary on a blank 8-bit, black and white, mask image.
		points = []
		pixels = []
		pts = geom.GetGeometryRef(0)
		
		for p in range(pts.GetPointCount()):
			points.append((pts.GetX(p), pts.GetY(p)))
		
		for p in points:
			pixels.append(self.world2Pixel(geoTrans, p[0], p[1])) # Transform nodes to geotrans of raster
		
		rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
		rasterize = ImageDraw.Draw(rasterPoly)
		rasterize.polygon(pixels, 0)
		
		# If the clipping features extend out-of-bounds and ABOVE the raster...
		if self.geoTrans[3] < maxY:
			# The clip features were "pushed down" to match the bounds of the
			#   raster; this step "pulls" them back up
			premask = self.imageToArray(rasterPoly)
			# We slice out the piece of our clip features that are "off the map"
			mask = numpy.ndarray((premask.shape[-2] - abs(iY), premask.shape[-1]), premask.dtype)
			mask[:] = premask[abs(iY):, :]
			mask.resize(premask.shape) # Then fill in from the bottom
			
			# Most importantly, push the clipped piece down
			geoTrans[3] = maxY - (maxY - gt[3])
		else:
			mask = self.imageToArray(rasterPoly)
		
		# Do the actual clipping
		if mask != None:
			try:
				clip2 = numpy.choose(mask,(clip, 0),mode='raise').astype(self.srcArray.dtype)
			except Exception:
				self.error = self.error + 1
				clip2 = None # Shape mismatch or Memory Error
			except ValueError:
				# Cut the clipping features to the raster
				rshp = list(mask.shape)
				if mask.shape[-2] != clip2.shape[-2]:
					rshp[0] = clip2.shape[-2]
				if mask.shape[-1] != clip2.shape[-1]:
					rshp[1] = clip2.shape[-1]
				# Resize to the clip
				mask.resize(*rshp, refcheck=False)

				try:
					clip2 = numpy.choose(mask,(clip,0),mode='raise').astype(self.srcArray.dtype)
				except ValueError:
					self.error = self.error + 1
					clip2 = None # Shape mismatch or Memory Error
		else:
			self.error = self.error + 1
			clip2 = None # Image to array failed because polygon outside range
		
		return clip2
	
	
	
	############################
	##### Active code
	
	
	### GET PARAMETERS
	
    #Define the folder address with the shapefile
	shpFolder = '/media/olivier/olivier_ext/TabulateAreaTEst' # args['pathshp']
	#Define the folder address with the image
	imgFolder = '/media/olivier/olivier_ext/TabulateAreaTEst' #args['pathtif']
	
	#Name of the shapefile
	shpName = 'grid_IA.shp' #args['shp']
	#Name of the image
	imgName = 'CDL2015_IA.tif' #args['tif']
	
	##No data value
	#if 'nodata' in args:
	#	no_data = args['nodata']
	#else:
	#	no_data = set()
	no_data = set(0)
	
	##Index to identify the polygons
	#if 'index' in args:
	#	index = args['index']
	#else:
	#	index = None
	index = 'GridCellID'
	
	
	### RUN ANALYSIS
	
	#Import the raster
	ds = gdal.Open(imgFolder+'/'+imgName)
	
	#Get the no data value and add with user-supplied if not the same
	if ds.GetRasterBand(1).GetNoDataValue():
		no_data.add(ds.GetRasterBand(1).GetNoDataValue())
	
	#Import the shapefile
	driver = ogr.GetDriverByName('ESRI Shapefile')
	dataSource = driver.Open(shpFolder+'/'+shpName, 0) # 0 means read-only. 1 means writeable.
	# Check to see if shapefile is found.
	if dataSource is None:
		print 'Could not open %s' % (shpFolder+'/'+shpName)
	else:
		layer = dataSource.GetLayer(0)
	
	#Crop the raster by the layer extent in case it is smaller
	args = {"INPUT":ds, "PROJWIN":", ".join(layer.GetExtent()), 
					"RTYPE":4, "COMPRESS":0, "BIGTIFF":3, "OUTPUT":None}
	if len(no_data) > 0:
		args["NO_DATA"] = no_data[0]
	ds = processing.runalg("qgis:cliprasterbyextent", args)
	ds = ds['OUTPUT'] #Get the address of the temporary file with the result of the algorithm
	ds = gdal.Open(ds) #Import again
	
	#Extract the raster projection information
	geo_transform = ds.GetGeoTransform()
	proj_ds = ds.GetProjection()
	
	#Get the data type for the array
	d_type = ds.dtype()
	
	#Convert the first band to an array
	try:
		ds_band = np.array(ds.GetRasterBand(1).ReadAsArray())
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
	rows, cols, n_bands = ds_band.shape()
	
	# Get Raster extent coordinates
	extent = GetExtentRaster(geo_transform,rows,cols)
	# Format to the same tuple structure as the vector layer
	extent = (extent[0][0],extent[2][0],extent[1][1],extent[0][1])
		
	
	#Prepare an id if the index is none
	i = 1
	
	tabulate = {}
	
	#Iterate over the features
	for feat in layer:
		
		# Get geometry of feature
		geom = feat.GetGeometryRef()
		#Loop through polygons if multipolygons
		if (geom.GetGeometryName() == 'MULTIPOLYGON'):
			U = []
			C = []
			count = 0
			for polygon in geom:
				geomInner = geom.GetGeometryRef(count)
				#Get polygon points
				ptsX, ptsY = getPtGeom(geomInner)   
				
				#Extract extent
				extent = [min(ptsX), max(ptsX), min(ptsY), max(ptsY)]
				
				#Crop the raster to the extent of the polygon
				clip = cropRasterArrayToExtent(ds_band, extent, geo_transform)
				
				if not clip:
					#This feature does not overlap with the raster
					continue
				
				#Transform the coordinates into pixels
				pixels = [world2Pixel(geo_transform, x, y) for x, y in zip(ptsX, ptsY)]
				
				#Get the pixel extent of the polygon
				pixExtent = world2Pixel(geo_transform, extent[0], extent[2]) +\
							world2Pixel(geo_transform, extent[1], extent[3])
				pixExtent = [pixExtent[0],pixExtent[2],pixExtent[1],pixExtent[3]]
				
				#Get the size of the rasterized polygon
				pxWidth = int(pixExtent[1] - pixExtent[0])
				pxHeight = int(pixExtent[3] - pixExtent[2])
				
				#Create a mask raster of the polygon
				rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
				rasterize = ImageDraw.Draw(rasterPoly)
				rasterize.polygon(pixels, 0)
				
				#Transform the image back to an array
				rasterPoly = imageToArray(rasterPoly)
				if not rasterPoly:
					return
				
				#Remove the pixels from this mask that do not overlap
				#with the raster
				#get the indices for trimming the array
				minX = min(0, -pixExtent[0])
				maxX = max(cols, pixExtent[1])-cols-1
				minY = max(rows, pixExtent[3])-rows
				maxY = -(min(0, -pixExtent[3]) + 1)
				#Trim
				rasterPoly = rasterPoly[minX:maxX,minY,maxY]
		
				#Apply the mask to the raster
				masked  = np.choose(rasterPoly,(clip, 0), mode='raise').astype(d_type)
				
				#Get unique values and counts
				u, c = np.unique(masked, return_counts=True)
				
				U.append(u)
				C.append(c)
				
				count += 1
			
			#Combine all the values of the individual polygons to put in tabulate
			u = set(x for l in U for x in l) #Get all unique values for the multipolygon
			#Loop through the combined values and get the sum of counts
			for uval in u:
				cval = 0
				for i in range(len(C)):
					try:
						ind = U[i].index(uval)
						cval += C[i][ind]
					except:
						continue
			
			#Export the combined values to tabulate
			
			
		elif (geom.GetGeometryName() == 'POLYGON'):
			#Get polygon points
			ptsX, ptsY = getPtGeom(geom)
			
			#Extract extent
			extent = [min(ptsX), max(ptsX), min(ptsY), max(ptsY)]
			
			#Crop the raster to the extent of the polygon
			clip = cropRasterArrayToExtent(ds_band, extent, geo_transform)
			
			if not clip:
				#This feature does not overlap with the raster
				continue
			
			#Transform the coordinates into pixels
			pixels = [world2Pixel(geo_transform, x, y) for x, y in zip(ptsX, ptsY)]
			
			#Get the pixel extent of the polygon
			pixExtent = world2Pixel(geo_transform, extent[0], extent[2]) +\
						world2Pixel(geo_transform, extent[1], extent[3])
			pixExtent = [pixExtent[0],pixExtent[2],pixExtent[1],pixExtent[3]]
			
			#Get the size of the rasterized polygon
			pxWidth = int(pixExtent[1] - pixExtent[0])
			pxHeight = int(pixExtent[3] - pixExtent[2])
			
			#Create a mask raster of the polygon
			rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
			rasterize = ImageDraw.Draw(rasterPoly)
			rasterize.polygon(pixels, 0)
			
			#Transform the image back to an array
			rasterPoly = imageToArray(rasterPoly)
			if not rasterPoly:
				return
			
			#Remove the pixels from this mask that do not overlap
			#with the raster
			#get the indices for trimming the array
			minX = min(0, -pixExtent[0])
			maxX = max(cols, pixExtent[1])-cols-1
			minY = max(rows, pixExtent[3])-rows
			maxY = -(min(0, -pixExtent[3]) + 1)
			#Trim
			rasterPoly = rasterPoly[minX:maxX,minY,maxY]
			
			#Apply the mask to the raster
			masked  = np.choose(rasterPoly,(clip, 0), mode='raise').astype(d_type)
			
			#Get unique values and counts
			u, c = np.unique(masked, return_counts=True)
			
			#Store unique values and counts in tabulate dictionary
			if not index:
				idf = i
			else:
				idf = feature.GetField(index)
            
			tabulate[idf] = (u, c)
            
            i += 1
            
	
	
	
	#Get the epsg code of the layer projection
	crs = layer.GetSpatialRef()
    crs.AutoIdentifyEPSG()
    code = srs.GetAuthorityCode(None)
    if code:
        srid= code
	
	
	for i in xrange(0,layer.GetFeatureCount()):
		poly = layer.GetFeature(i)
		#if self.lyr.GetFeatureCount() == 1:
		#	poly = layer.GetFeature(0)
		#else:
		#	poly = layer.GetFeature(i)

		# Test if polygon feature is inside raster extent, otherwise return None as result
		array = getClipArray(poly)
		if array is None: # Multi polygon or no raster values below ?
			next
                
		classes = sorted(numpy.unique(array)) # get classes
		for val in (self.nodata,0):# Remove raster nodata value and zeros from class list
			try:
				classes.remove(val)
			except ValueError:
				pass
	
	
	
	
	#Rasterize the shapefile
	driver = gdal.GetDriverByName('MEM')  # In memory dataset
    mask = driver.Create('', cols, rows, 1, gdal.GDT_UInt16)
    mask.SetGeoTransform(geo_transform)
    mask.SetProjection(projection)
    gdal.RasterizeLayer(mask, [1], layer, options = ["ATTRIBUTE="+index])
    
	mask_band = np.array(mask.GetRasterBand(1).ReadAsArray())
	mask = None #Close the raster to free memory space
	
	#Combine the two arrays into one
	
	
	#Loop through the features
	
	
	
	#Remove any layers with same name in registry if already one
	for layer in QgsMapLayerRegistry.instance().mapLayers().values():
		if layer.name() == "Polygons" or layer.name() == "Image":
			QgsMapLayerRegistry.instance().removeMapLayers([layer.id()])
					
	#Import the ground shapefile of polygons
	vt = iface.addVectorLayer(shpFolder+"/"+shpName, "Polygons", "ogr")
	if not vt or not vt.isValid():
			print("Failed to load %s" %(shpFolder+"/"+shpName))
			sys.exit()
            
	#Import the raster image
	img = iface.addRasterLayer(imgFolder+"/"+imgName, "Image")
	if not img or not img.isValid():
			print("Failed to load %s" %(imgFolder+"/"+imgName))
			sys.exit()
	
	#Create an empty list to hold the index of the polygons, 
	#Unique values and counts
	tabulate = []
	
	#Add an index to be able to identify the vt polygons uniquely
	#if no index was specified by the user.
	i = 1
	
    #Loop over all polygons in the VT
	for f in vt.getFeatures():
		
		#Select the current polygon
		vt.select(f.id())
		#Export to a temporary file
		QgsVectorFileWriter.writeAsVectorFormat(vt, shpFolder+"/temp.shp", "utf-8", 
						None, "ESRI Shapefile", True) #Last 1 boolean is for exporting selection only
		"""Can I improve speed by using layer in memory? Would mean the processing needs to accept in memory"""
		
		#Add it to the display window to check
		#QgsMapLayerRegistry.instance().addMapLayer(F)
		
		#Clip the raster by the current polygon
		if not no_data is None:
			first = processing.runalg("gdalogr:cliprasterbymasklayer", 
					{"INPUT":img, "MASK":shpFolder+"/temp.shp",  
					"NO_DATA":no_data, "CROP_TO_CUTLINE":1, "KEEP_RESOLUTION":1,
					"COMPRESS":0, "RTYPE":5, "BIGTIFF": 2,
					"TFW":0, "OUTPUT":None})
		else:
			first = processing.runalg("gdalogr:cliprasterbymasklayer", 
					{"INPUT":img, "MASK":shpFolder+"/temp.shp",  
					"CROP_TO_CUTLINE":1, "KEEP_RESOLUTION":1,
					"COMPRESS":0, "RTYPE":5, "BIGTIFF": 2,
					"TFW":0, "OUTPUT":None})
		first = first['OUTPUT']
		first = QgsRasterLayer(first, "clip")
		
		#Extract all the non null values from the raster
		prov = first.dataProvider()
		width = prov.xSize()
		height = prov.ySize()
		block = prov.block(1, prov.extent(), width, height)
		
		values = []
		for x in range(width):
			for y in range(height):
				values.append(int(block.value(x, y)))
		
		#Remove no data values if specified
		#Pass the data to integer
		if no_data is not None:
			values = [str(v) for v in values if not v == no_data]
		
		#Get the unique values
		values = Counter(values)
		
		#Add the id of polygon
		if index == None:
			values['id'] = f.id()
		else:
			values['id'] = f[index]
		
		#Add the results to the list of tabulated values
		tabulate.append(values)
		
		#Update the index
		i += 1
		#Remove selection before next iteration
		vt.removeSelection()
		
	#Remove temporary files
	if os.path.exists(shpFolder+"/temp.shp"):
		ogr.GetDriverByName("ESRI Shapefile").DeleteDataSource(shpFolder+"/temp.shp")
	if os.path.exists(imgFolder+"/temp.tif"):
		os.remove(imgFolder+"/temp.tif")
	if os.path.exists(imgFolder+"/temp.tif.aux.xml"):
		os.remove(imgFolder+"/temp.tif.aux.xml")
	
	
	print(' '.join('{}--{}'.format(key, val) for key, val in tabulate[0].items()))
	print(' ')
	print(' '.join(tabulate[0].keys()))
	print(' ')
	print(' '.join(tabulate[28].keys()))
	#Add the missing keys for all polygons
	allkeys = []
	for d in tabulate:
		allkeys = allkeys + d.keys()
	allkeys = set(allkeys)
	#set().union(*(d.keys() for d in tabulate))
	
	print(' '.join(sorted(allkeys)))
	for i in tabulate:
		for missing in allkeys.difference(i):
			i[missing] = 0
	
	#Export list of dictionaries to a tab delimited file
	with open(shpFolder+'/test.txt', "w") as f:
		dict_writer = DictWriter(f, ['id']+sorted(allkeys).remove('id'), delimiter="\t")
		dict_writer.writeheader()
		for p in tabulate:
			dict_writer.writerow(p)
