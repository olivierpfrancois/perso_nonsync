# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

"""
The algorithm imports the raster and transforms it into a numpy 
array. It loops through the polygons in the vector layer, crops the 
raster to the polygon, transforms the polygon into a numpy array (using 
Image and ImageDraw) using the pixel coordinates of the nodes of the 
polygon.
Then it uses the polygon array to mask the raster and counts the unique
values of the pixels inside the polygon.

The tabulated values are exported into a tab delimited file.
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
import Image, ImageDraw
import shapely
from shapely import wkt


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
		if len(srcArray.shape) == 2:
			rows, cols = srcArray.shape
		else:
			rows, cols, n_bands = srcArray.shape
		
		#If all positions are out of the array, return None
		if max(lrX,lrY) == 0 or (ulY >= rows and ulX >= cols):
			return 

		# Clip the raster to the shapes bounding box
		clip = srcArray[ulY:(lrY+1), ulX:(lrX+1)]
		
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
		except:
			try:
				a = np.fromstring(i.tobytes(),'b')   
			except SystemError:
				return None
		if a != None:
			a.shape=i.im.size[1], i.im.size[0]
		return a
	
	
	############################
	##### Active code
	
	
	### GET PARAMETERS
	
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
	
	#Get the data type for the array
	d_type = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
	if not d_type in ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]:
		print("The data type of the raster is not supported by the algorithm")
		return
	
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
		
		d_types = ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]
		rtype = d_types.index(d_type)
		
		#Crop the raster by the layer extent in case it is smaller
		args = {"INPUT":imgFolder+'/'+imgName, "PROJWIN":", ".join(map(str, layer.GetExtent())), 
					"RTYPE":rtype, "COMPRESS":0, "BIGTIFF":3, "OUTPUT":None}
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
		
		#Loop through polygons if multipolygons
		if (geom.GetGeometryName() == 'MULTIPOLYGON'):
			U = []
			C = []
			count = 0
			for polygon in geom:
				geomInner = geom.GetGeometryRef(count)
				if geomInner.GetPointCount() == 0:
					geomInner = geomInner.GetGeometryRef(0)
				
				
				#Get unique values and counts
				#u, c = np.unique(masked, return_counts=True)
				
				U.append(u)
				C.append(c)
			
				count += 1
			
			#Combine all the values of the individual polygons to put in tabulate
			u = list(set(x for l in U for x in l)) #Get all unique values for the multipolygon
			c = []
			#Loop through the combined values and get the sum of counts
			for uval in u:
				cval = 0
				for i in range(len(C)):
					try:
						ind = U[i].index(uval)
						cval += C[i][ind]
					except:
						continue
				c.append(cval)
			
			#Export the combined values to tabulate
			#Store unique values and counts in tabulate dictionary
            
			tabulate[idf] = (u, c)
			uniques.update(u)
			
		elif (geom.GetGeometryName() == 'POLYGON'):
			
			#Transform polygon into shapely object
			poly = wkt.loads(geom.ExportToWkt())
			
			if geom.GetPointCount() == 0:
				geom = geom.GetGeometryRef(0)
			
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
			
			if ((pixExtent[0] > rows and pixExtent[2] > cols) or 
				(pixExtent[1] == 0 and pixExtent[3] == 0)):
				#Polygon is outside of raster
				continue
			
			#crop the raster based on this extent
			clip = ds_band[pixExtent[2]:pixExtent[3],pixExtent[0]:pixExtent[1]]
			
			#Get polygon points
			pts = []
			for i in range(geom.GetPointCount()):
				# GetPoint returns a tuple
				pt = geom.GetPoint(i)
				#Transform into pixel and line positions
				pts.append(world2Pixel(geo_transform, pt[0], pt[1]))
			#Unpack the pixel positions
			ptsX, ptsY = map(list,zip(*pts))
			
			#Get the size of the rasterized polygon
			pWidth = int(max(ptsX)+1 - min(ptsX))
			pHeight = int(max(ptsY)+1 - min(ptsY))
			
			#Create a mask raster of the polygon
			pixels = [(t[0]-min(ptsX),t[1]-min(ptsY)) for t in pts] #Reset pixel values
			#print("pixels: "+' '.join(['_'.join(map(str,x)) for x in pixels]))
			rasterPoly = Image.new("L", (pWidth, pHeight), 1) #Create image of 1
			rasterize = ImageDraw.Draw(rasterPoly) #Prepare drawing object
			#Create polygon of 0s inside and 2 on the outline
			rasterize.polygon(pixels, fill=0, outline=2) 
			
			#Transform the image back to an array
			rasterPoly = imageToArray(rasterPoly)
			if rasterPoly == None:
				return
			
			#Get the pixels that have a value of 2 to test if inside polygon
			for x, y in np.ndindex(rasterPoly.shape):
				if rasterPoly[x,y] == 2:
					#Reset position as in original raster
					L = x + min(ptsY)
					P = y + min(ptsX)
					
					#Transform into coordinates of the center of the pixel
					Xp = xOrigin + P*pixelWidth + L*xRotation + pixelWidth/2
					Yp = yOrigin + P*yRotation + L*pixelHeight + pixelHeight/2
					
					#Transform into shapely point
					pt = shapely.geometry.Point(Xp,Yp)
					
					#Check if point inside the polygon
					if pt.within(poly):
						rasterPoly[x,y] = 0 #Set it to 0 to consisder it
					else:
						rasterPoly[x,y] = 1 #Set it to 1 to discard it
			
			#Remove the pixels from this mask that do not overlap
			#with the raster
			#get the indices for trimming the array
			minP = 0 + min(0, min(ptsX)) #If pixel value is >0, inside raster
			maxP = pWidth - (max(cols, max(ptsX)+1)-cols)
			minL = 0 + min(0, min(ptsY))
			maxL = pHeight - (max(rows, max(ptsY)+1)-rows)
			#Trim
			rasterPoly = rasterPoly[minL:maxL,minP:maxP]
			
			#Apply the mask to the raster
			#Where raster mask == 0 take values from clip, where mask == 1 take 0
			if len(no_data) > 0:
				holder = no_data[0]
			else:
				holder = 0
			masked  = np.choose(rasterPoly,(clip, holder), mode='raise')
			
			#Get unique values and counts
			u, c = np.unique(masked, return_counts=True)
			#Remove no data
			for n in no_data:
				c = c[u!=n]
				u = u[u!=n]
			
			#Store unique values and counts in tabulate dictionary
			tabulate[id_in] = dict(zip([prefix+x for x in map(str,u)], c))
	
	#Transform the result dictionary for the export
	out = []
	#Prepare a set to hold the unique values in the raster
	uniques = set()
	for k, v in tabulate.items():
		uniques.update(v.keys())
		v['id'] = k
		out.append(v)
	
	#Export to a tab delimited file
	with open(shpFolder+'/test.txt', "w") as f:
		dict_writer = DictWriter(f, ['id']+sorted(uniques, key = lambda x: float(x.split("_")[1])), delimiter="\t", restval="0")
		dict_writer.writeheader()
		for p in out:
			dict_writer.writerow(p)
	
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
