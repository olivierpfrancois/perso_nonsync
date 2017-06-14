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
import qgis.utils

from osgeo import ogr, gdal
from csv import DictWriter
import os, processing, re
import numpy as np
import time

class Tabulate:
	def __init__(self, iface):
		"""Initialize using the qgis.utils.iface
		object passed from the console.
		"""
		self.iface = qgis.utils.iface
		
	def GetExtentRaster(self, gt, cols, rows):
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
	
	def new_raster_from_base(self, base, outputURI, format, nodata, datatype):
		''' Create an empty copy of a raster from an existing one
			@type base:   gdal raster layer
			@param base: raster
			@type outputURI:   string
			@param outputURI: address + name of the output raster (extension should agree with format, none for memory)
			@type format:   string
			@param format: format for the dataset (e.g. "GTiff", "MEM")
			@type nodata:   int/float
			@param nodata: no data value (type should agree with raster type)
			@type datatype:   gdal data type (e.g. gdal.GDT_Int32)
			@param datatype: data type for the raster
			@rtype:    gdal raster layer
			@return:   raster filled with no data values
		'''
		
		cols = base.RasterXSize
		rows = base.RasterYSize
		projection = base.GetProjection()
		geotransform = base.GetGeoTransform()
		bands = base.RasterCount

		driver = gdal.GetDriverByName(format)

		new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype)
		new_raster.SetProjection(projection)
		new_raster.SetGeoTransform(geotransform)

		for i in range(bands):
			new_raster.GetRasterBand(i + 1).SetNoDataValue(nodata)
			new_raster.GetRasterBand(i + 1).Fill(nodata)

		return new_raster
	
	def tabulate_area(self, rasterName, bands, shpName, index=None, tabulate=True, prefix='', statistics=None, outTxt=None, alltouch=False):
		'''
		Computes the tabulate area for each polygon in the shapefile on each of the required bands
		
		INPUTS:
		rasterName: string
			Full address of the raster including extension
			
		bands: list of int
			list of band numbers to extract information from
		
		shpName: string
			Full address of the shapefile including extension
		
		index: [optional] string
			Name of the variable in the shapefile to use as labels for the polygons
		
		tabulate: [optional] logical
			If True the function will do a tabulate for each of the polygons
		
		prefix: [optional] string
			prefix to use to write in front of the raster values in the output file and dictionaries
		
		statistics: [optional] list of string
			The function will return for each polygon the statistics listed in this option. 
			The argument should be None if no statistics are wanted
			The statistics can be chosen among:
			'min','max','mean','count','sum','sd','median','unique'
			Where count is the number of pixels in the polygon, and 
			unique the number of unique values in the polygon
			THE STATISTICS ARE ONLY EXPORTED TO A TXT FILE.
		
		outTxt: [optional] string
			Full address of the output .txt including extension
			
		alltouch: [optional] logical
			Parameter value for the rasterization.
			If True, all pixels touched by the polygons will be considered as part of the polygons.
			If False, only those with the centroid inside the polygons will be considered.
			Defaults to False
			
		The function returns a list, one element for each band.
		Each element is a dictionary of dictionaries has the label of the polygon.
		The main dictionary has one element per polygon, where the polygon label is the key.
		The sub-dictionaries have one element per value found in the raster band, 
		where the value is the key and the count is the value.
		'''
		
		if not tabulate and not statistics:
			print("No tabulation or statistics specified.")
			return
		
		#Transform alltouch
		if alltouch:
			alltouch = 'TRUE'
		else:
			alltouch = 'FALSE'
		
		#Import the raster
		try:
			ds = gdal.Open(rasterName)
		except RuntimeError, e:
			print 'Unable to open '+rasterName
			print e
			return
		
		#Get the number of bands
		nBands = ds.RasterCount
		
		#Get the no data value
		no_data = ds.GetRasterBand(1).GetNoDataValue()
		
		#Get the data type for the array
		d_type = gdal.GetDataTypeName(ds.GetRasterBand(1).DataType)
		if not d_type in ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]:
			print("The data type of the raster is not supported by the algorithm")
			return
		
		#Import the shapefile
		driver = ogr.GetDriverByName('ESRI Shapefile')
		dataSource = driver.Open(shpName, 1) # 0 means read-only. 1 means writeable.
		# Check to see if shapefile is found.
		if dataSource is None:
			print 'Could not open %s' % (shpName)
			return
		else:
			layer = dataSource.GetLayer(0)
		
		#Add an id field for the rasterization
		layerNames = [layer.GetLayerDefn().GetFieldDefn(n).GetName() for n in range(layer.GetLayerDefn().GetFieldCount())]
		if not "IDRAST" in layerNames:
			idField = ogr.FieldDefn("IDRAST", ogr.OFTInteger)
			layer.CreateField(idField)
		idField = layer.GetLayerDefn().GetFieldIndex("IDRAST")
		
		#Get the index of the field to label values
		if index:
			#Get the index of the field with the polygon ids
			idx = layer.GetLayerDefn().GetFieldIndex(index)
		else:
			idx = idField
		
		#Prepare the labels for the output file
		labels = []
		i = 1
		for feat in layer:
			feat.SetField("IDRAST", i)
			layer.SetFeature(feat)
			i += 1
			labels.append(feat.GetFieldAsString(idx))
		ids = range(1,i)
		
		#Get the extent of the shapefile layer
		shp_ext = layer.GetExtent()
		
		# Get Raster extent coordinates
		rast_ext = self.GetExtentRaster(ds.GetGeoTransform(), ds.RasterYSize, ds.RasterXSize)
		# Format to the same tuple structure as the vector layer
		rast_ext = (rast_ext[0][0],rast_ext[2][0],rast_ext[1][1],rast_ext[0][1])
		'''
		#Compare the two extent to see if raster inside the shapefile
		if (rast_ext[0] < shp_ext[0] or rast_ext[1] > shp_ext[1] or 
			rast_ext[2] < shp_ext[2] or rast_ext[3] > shp_ext[3]):
			
			d_types = ["Byte","Int16","UInt16","UInt32","Int32","Float32","Float64"]
			rtype = d_types.index(d_type)
			
			#Crop the raster by the layer extent in case it is smaller
			args = {"INPUT":rasterName, "PROJWIN":", ".join(map(str, layer.GetExtent())), 
						"RTYPE":rtype, "COMPRESS":0, "BIGTIFF":3, "OUTPUT":None}
			if no_data:
				args["NO_DATA"] = no_data
			
			ds = processing.runalg("gdalogr:cliprasterbyextent", args)
			ds = ds['OUTPUT'] #Get the address of the temporary file with the result of the algorithm
			ds = gdal.Open(ds) #Import again
		'''
		#Prepare an empty raster to rasterize the shapefile
		#rVector = self.new_raster_from_base(ds, 'temp.tif', 'GTiff', -1, gdal.GDT_Int32)
		rVector = self.new_raster_from_base(ds, 'temp', 'MEM', -1, gdal.GDT_Int32)
		
		#Rasterize the vector layer
		gdal.RasterizeLayer(rVector, [1], layer, options=['ALL_TOUCHED='+alltouch, 'ATTRIBUTE=IDRAST'])
		
		#Remove rasterization id
		layer.DeleteField(idField)
		
		#Remove layer
		del layer
		del dataSource
		
		try:
			v_band = np.array(rVector.GetRasterBand(1).ReadAsArray())
			if len(v_band.shape) == 3:
				v_band = v_band[:,:,0] #Keep only the first two dimensions
			
		except ValueError:
			print("Raster file is too big for processing. Please crop the file and try again.")
			
			#Close the rasters
			ds = None
			rVector = None
			
			return
		
		#Close the raster
		rVector = None
		
		#Flatten the array
		v_band.ravel()
		
		#Get no data positions
		remove_nodata = v_band!=-1
		
		#Remove the no data
		v_band = v_band[remove_nodata]
		
		#Order the data
		order = np.argsort(v_band)
		v_band = v_band[order]
		
		#Get the indices of the polygon ids
		pids = [np.searchsorted(v_band, ind) for ind in ids]
		pids.append(-1)
		
		#Prepare an empty list to hold the dictionaries for each band
		tabout = []
		
		#Prepare an empty list to hold the dictionaries for each band
		stats = []
				
		if statistics:
			statistics = [s for s in statistics if s in ['min','max','mean','count','sum','sd','median','unique']]
			statistics.sort()
			
		for b in bands:
			#Convert the band to an array
			try:
				ds_band = np.array(ds.GetRasterBand(b).ReadAsArray())
				if len(ds_band.shape) == 3:
					ds_band = ds_band[:,:,0] #Keep only the first two dimensions
				
				#Flatten the array
				ds_band.ravel()
				
			except ValueError:
				print("Raster file is too big for processing. Please crop the file and try again.")
				
				#Close the rasters
				ds = None
				
				return
			
			#Remove no data values and other values
			ds_band = ds_band[remove_nodata]
			
			#Order the data
			ds_band = ds_band[order]
			
			if tabulate:
				#Prepare a dictionary to hold the counts
				tabout.append({})
				
				for ind in range(len(ids)):
					
					#Get unique values and counts
					u, c = np.unique(ds_band[pids[ind]:pids[ind+1]], return_counts=True)
					
					#Store unique values and counts in tabulate dictionary
					tabout[-1][labels[ind]] = dict(zip([prefix+x for x in map(str,u)], c))
			
			if statistics:
				#Prepare a dictionary to hold the counts
				stats.append({})
				
				for ind in range(len(ids)):
					
					values = ds_band[pids[ind]:pids[ind+1]]
					
					res = {}
					
					#Get the statistics
					if 'min' in statistics:
						res['min'] = np.amin(values)
					if 'max' in statistics:
						res['max'] = np.amax(values)
					if 'mean' in statistics:
						res['mean'] = np.mean(values)
					if 'count' in statistics:
						res['count'] = values.size
					if 'sum' in statistics:
						res['sum'] = np.sum(values)
					if 'sd' in statistics:
						res['sd'] = np.std(values)
					if 'median' in statistics:
						res['median'] = np.median(values)
					if 'unique' in statistics:
						u = np.unique(values)
						res['unique'] = u.size
					
					#Store statistics in stats dictionary
					stats[-1][labels[ind]] = res
		
		if outTxt:
			if tabulate:
				for b, d in zip(bands, tabout):
					#Transform the result dictionary for the export
					out = []
					#Prepare a set to hold the unique values in the raster
					uniques = set()
					for k, v in d.items():
						uniques.update(v.keys())
						v['id'] = k
						out.append(v)
					
					#Remove the no data value if any
					if no_data:
						uniques.discard(n) 
					
					#Export to a tab delimited file
					with open(re.sub(".txt","_"+str(b)+".txt",outTxt), "w") as f:
						dict_writer = DictWriter(f, ['id']+sorted(uniques, key = lambda x: float(x.split("_")[1])), extrasaction='ignore', delimiter="\t", restval="0")
						dict_writer.writeheader()
						for p in out:
							dict_writer.writerow(p)
			
			if statistics:
				for b, s in zip(bands, stats):
					#Transform the result dictionary for the export
					out = []
					#Prepare a set to hold the unique values in the raster
					for k, v in s.items():
						v['id'] = k
						out.append(v)
					
					#Export to a tab delimited file
					with open(re.sub(".txt","_STATS_"+str(b)+".txt",outTxt), "w") as f:
						dict_writer = DictWriter(f, ['id']+statistics, extrasaction='ignore', delimiter="\t", restval="0")
						dict_writer.writeheader()
						for p in out:
							dict_writer.writerow(p)
			
			
		#Close the rasters
		ds = None
			
		#Return the list of dictionaries
		if tabulate:
			return tabout
	
def run_script(iface):
	
	t0 = time.time()
	
	tabu = Tabulate(iface)
	
    ### GET PARAMETERS
    #Define the folder address with the shapefile
	shpFolder = 'D:/gedata_current/jde_coffee/MODIS/ES' # args['pathshp'] #'/media/olivier/olivier_ext/TabulateAreaTEst'
	#Define the folder address with the image
	imgFolder =  'D:/gedata_current/jde_coffee/MODIS/ES/smooth_data' #args['pathtif'] #'/media/olivier/olivier_ext/TabulateAreaTEst'
	
	#Name of the shapefile
	shpName = 'ES_areas_segments_modis_sum.shp' #'smallpols.shp' #args['shp']
	#Name of the image
	imgName = 'smooth_MOD13Q1_2005-01-01_h13-14v10-11.250m_16_days_NDVI.tif' #args['tif']
	
	##Index to identify the polygons
	#if 'index' in args:
	#	index = args['index']
	#else:
	#	index = None
	index = 'serp_id' #'ID' 
	
	##Prefix for the data
	prefix = 'cat_'
	
	##Option for the rasterization
	alltouch = 'FALSE'
	
	### RUN ANALYSIS
	tabu.tabulate_area(rasterName=imgFolder+'/'+imgName, bands=[1], shpName=shpFolder+'/'+shpName, 
							index=index, tabulate=True, prefix=prefix, statistics=[],
							outTxt=shpFolder+'/test.txt', alltouch=alltouch)
	
	t1 = time.time()
	print(str(t1-t0))