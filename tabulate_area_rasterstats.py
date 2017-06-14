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

import numpy as np
import warnings, os
from affine import Affine
from shapely.geometry import shape
from rasterstats.io import read_features, Raster
from rasterstats.utils import (rasterize_geom, get_percentile,
                    remap_categories, key_assoc_val, boxify_points)
from csv import DictWriter
import time

def run_script(iface):

	def gen_tabulate(
		vectors, 
		raster,
		layer=0,
		index=None,
		band_num=1,
		nodata=None,
		affine=None,
		all_touched=False,
		categorical=False,
		category_map=None,
		prefix=None):
		"""Zonal statistics of raster values aggregated to vector geometries.
		
		Parameters
		----------
		vectors: path to an vector source or geo-like python objects
		
		raster: ndarray or path to a GDAL raster source
			If ndarray is passed, the ``affine`` kwarg is required.
		
		layer: int or string, optional
			If `vectors` is a path to an fiona source,
			specify the vector layer to use either by name or number.
			defaults to 0
		
		index: string, optional
			The name of the variable in the vector shapefile that will be
			used to id the polygons in the output file
		
		band_num: int, optional
			If `raster` is a GDAL source, the band number to use (counting from 1).
			defaults to 1.
		
		nodata: float, optional
			If `raster` is a GDAL source, this value overrides any NODATA value
			specified in the file's metadata.
			If `None`, the file's metadata's NODATA value (if any) will be used.
			defaults to `None`.
		
		affine: Affine instance
			required only for ndarrays, otherwise it is read from src
		
		all_touched: bool, optional
			Whether to include every raster cell touched by a geometry, or only
			those having a center point within the polygon.
			defaults to `False`
		
		categorical: bool, optional
		
		category_map: dict
			A dictionary mapping raster values to human-readable categorical names.
			Only applies when categorical is True
		
		prefix: string
			add a prefix to the keys (default: None)
		
		Returns
		-------
		generator of dicts
			Each item corresponds to a single vector feature and
			contains keys for each of the specified stats.
		"""

		# Handle 1.0 deprecations
		transform = kwargs.get('transform')
		if transform:
			warnings.warn("GDAL-style transforms will disappear in 1.0. "
						  "Use affine=Affine.from_gdal(*transform) instead",
						  DeprecationWarning)
			if not affine:
				affine = Affine.from_gdal(*transform)
		
		tabulate = {}
		
		with Raster(raster, affine, nodata, band_num) as rast:
			features_iter = read_features(vectors, layer)
			for i, feat in enumerate(features_iter):
				#Get the index for the geometry
				if index:
					id_in = feat['properties'][index]
				else:
					id_in = feat['id']
				
				geom = shape(feat['geometry'])

				if 'Point' in geom.type:
					geom = boxify_points(geom, rast)

				geom_bounds = tuple(geom.bounds)

				fsrc = rast.read(bounds=geom_bounds)

				# create ndarray of rasterized geometry
				rv_array = rasterize_geom(geom, like=fsrc, all_touched=all_touched)
				assert rv_array.shape == fsrc.shape

				# Mask the source data array with our current feature
				# we take the logical_not to flip 0<->1 for the correct mask effect
				# we also mask out nodata values explicitly
				masked = np.ma.MaskedArray(
					fsrc.array,
					mask=np.logical_or(
						fsrc.array == fsrc.nodata,
						np.logical_not(rv_array)))

				if masked.compressed().size == 0:
					# nothing here, fill with zeros and move on
					continue
				else:
					keys, counts = np.unique(masked.compressed(), return_counts=True)
					if prefix:
						key = [prefix+str(k) for k in keys]
					#Add the counts to the tabulated values as a dictionnary
					tabulate[id_in] = dict(zip(keys, counts))
				
		return(tabulate)

	
	def export_tabulate(tabulate, outName):
		#Transform the result dictionary for the export
		out = []
		#Prepare a set to hold the unique values in the raster
		uniques = set()
		for k, v in tabulate.items():
			uniques.update(v.keys())
			v['id'] = k
			out.append(v)
		
		#Export to a tab delimited file
		with open(outName, "w") as f:
			dict_writer = DictWriter(f, ['id']+sorted(uniques, key = lambda x: int(x.split("_")[1])), 
					delimiter="\t", restval="0")
			dict_writer.writeheader()
			for p in out:
				dict_writer.writerow(p)
	
	
	
	#Define the folder address with the shapefile
	shpFolder = 'D:\TabulateAreaTEst' #'/media/olivier/olivier_ext/TabulateAreaTEst' # args['pathshp']
	#Define the folder address with the image
	imgFolder = 'D:\TabulateAreaTEst' #'/media/olivier/olivier_ext/TabulateAreaTEst' #args['pathtif']
	
	#Name of the shapefile
	shpName = 'grid_IA.shp' #args['shp']
	#Name of the image
	imgName = 'CDL2015_IA.tif' #args['tif']
	
	t0 = time.time()
	
	tabulate = gen_tabulate(
		vectors=os.path.join(shpFolder,shpName), 
		raster=os.path.join(imgFolder,imgName),
		layer=0,
		band_num=1,
		nodata=None,
		affine=None,
		all_touched=False,
		categorical=False,
		category_map=None,
		prefix='cat_')
	
	export_tabulate(tabulate, os.path.join(shpFolder,"test_rasterio.txt"))
    
	t1 = time.time()
	print(str(t1-t0))
	
