# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" Rasterizes a shapefile to a wanted resolution

    The script creates an empty raster at the wanted resolution 
    in the bounding box of the shapefile (using the nodata value provided), 
    and then transfers the values of the selected attribute based on where the 
    pixel centroid falls compared to the shapefile polygons. 
    
    Variables:
    shp: Full address of the polygon shapefile to rasterize
    attrib: Name of the attribute to use in the shapefile to populate the raster pixels
    resolution: wanted resolution for the raster, in meters if UTM
    nodata: No data value to use for the raster
    
    The output tiff is in the same folder as the shapefile, with the same name.
    
     """

# Some commonly used imports

import os, re
from osgeo import gdal, ogr
    
def toRaster(nm, attrib, resRaster, nodata):
    
    # Open the data source and read in the extent
    source_ds = ogr.Open(nm)
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
    #Create the name of the output raster
    outName = re.sub(".shp",".tif",nm)
    
    #Remove output file if it already exists
    try:
        os.remove(outName)
    except OSError:
        pass
    
    # Create the destination data source
    x_res = int((x_max - x_min) / resRaster)
    y_res = int((y_max - y_min) / resRaster)
    target_ds = gdal.GetDriverByName('GTiff').Create(outName, x_res, y_res, 1, gdal.GDT_UInt16, options=['COMPRESS=LZW'])
    target_ds.SetProjection(source_layer.GetSpatialRef().ExportToWkt())
    target_ds.SetGeoTransform((x_min, resRaster, 0, y_max, 0, -resRaster))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, options=['ALL_TOUCHED=FALSE', 'ATTRIBUTE='+attrib])
    
    target_ds = None
    source_ds = None


if __name__ == '__main__':
    
    #PARAMETERS
    
    #Address of the shapefile to rasterize
    shp = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/Vietnam/LD_grid_c_m/LD_grid_c_m.shp'
    
    #Output resolution wanted for the raster (in meters if utm)
    resolution = 6

    #Attribute to use for the raster values
    attrib = 'ID'
    
    #Do the rasterization
    toRaster(shp, attrib, resolution)