from osgeo import ogr, gdal, osr, gdalconst
from csv import DictWriter
import os, processing, re, csv, platform, sys
import numpy as np

sys.path.append('C:/Users/Olivier/OlivierGithub/QGIS-scripts')
import image_proc_functions as img

def run_script(iface):
    ##########################################################################################################################
    ##########################################################################################################################
    ##PARAMETERS
    
    #Address of the folder with the classifications
    classifFolder = 'E:/gedata_current/jde_coffee/data/SDM/classifs'
    #TAB delimited file with the priority order in which to combine the classifications
    #The file should have two columns, the first is a priority number, the second the file name of the classification
    orderClassif = 'E:/gedata_current/jde_coffee/data/SDM/classifs/classif_priorities_SDM.txt'
    #Classes to pass to NA values in combining the classifications
    #Leave an empty list if None
    #These names should not include any trailing numbers if several categories have been made (e.g. Cloud1, Cloud2...)
    toNa = ['Cloud','Shadow']
    #Additional values in the classifications that are not in the legends
    valuesNotInLegend = [0]
    
    
    #Address of the folder with the legends
    legendFolder = 'E:/gedata_current/jde_coffee/data/SDM/classifs/Legend'
    #Legend file names prefix to add to the classification names
    legendPrefix = 'legend_'
    #Legend file extensions: could be '.csv' or '.txt'
    legendExt = '.csv'
    #Legend file delimiter
    legendDelim = ','
    
    #Address of the folder with the confusion matrices
    matrixFolder = 'E:/gedata_current/jde_coffee/data/SDM/classifs/Matconf'
    #Matrix file names prefix to add to the classification names
    matrixPrefix = 'Matconf_'
    #Matrix file extensions: could be '.csv' or '.txt'
    matrixExt = '.csv'
    #Legend file delimiter
    matrixDelim = ','
    
    #Address of the folder for the temporary files
    tempFolder = 'E:/gedata_current/temp'
    
    
    #Address and name of the grid shapefile
    gridShape = 'E:/gedata_current/jde_coffee/data/SDM/grids/AOI_grid_SDM_2km.shp'
    #Field to use for the ID of the cells
    idField = 'serp_id'
    
    
    #Projection epsg for the output. 
    #Leave None to use the one from the grid. Should be the same as the grid if the grid doesn't have a projection info!
    epsgProject = 32723
    
    
    #SRTM Mask option
    srtmMask = False
    srtmRaster = 'E:/gedata_current/jde_coffee/data/SDM/srtm/srtm_v3_SDMnet.tif'
    breakValues = [550] #list of values used as break to separate the classifications 
    labelsElev = ['robusta','arabica']
    
    
    #Share of clouds/missing per cell under which correct by simple proportional scaling
    propCloud <- 0.5
    
    
    #Output name and address for the combined files and legend
    exportClassifName = 'megot_classif.tif'
    
    
    
    ##########################################################################################################################
    ##########################################################################################################################
    ##CODE
    
    #Combine the classifications without and with the clouds
    combine_classifications(classifFolder, orderClassif, legendFolder, legendPrefix, legendExt, 
                            legendDelim, tempFolder, 'nocloud_'+exportClassifName, toNa=toNa, valuesNotInLegend=valuesNotInLegend)
    
    combine_classifications(classifFolder, orderClassif, legendFolder, legendPrefix, legendExt, 
                            legendDelim, tempFolder, 'cloud_'+exportClassifName, toNa=[], valuesNotInLegend=valuesNotInLegend)
    
    #Merge the two combined rasters
    #Check the computer platform
    ordi = checkPlatform()
    if ordi == "Windows":
        import gdal_merge as gm
    #Prepare the arguments for gdal_merge.py -- They all need to be strings
    args = ['-n', '255', '-a_nodata', '255', '-of', 'GTiff', '-ot', 'UInt32', '-co', 'COMPRESS:LZW', '-o', classifFolder+'/'+exportClassifName]
    args.extend(['cloud_'+exportClassifName,'nocloud_'+exportClassifName])
    #Do the merge
    if ordi == "Linux":
        args = ['gdal_merge.py','-of', "GTiff"]+args
        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        p.communicate()[0] #Waits until the subprocess is finished before running the rest of the code
    elif ordi == "Windows":
        #Windows makes an error when calling gdal_merge directly, 
        #So need to specify python for the executable and to make it run in the shell
        sys.argv = args
        gm.main()
    
    #Combine the legends to have a single legend?
    
    #Remove the intermediary files
    for f in [classifFolder+'/cloud_'+exportClassifName,classifFolder+'/nocloud_'+exportClassifName,legendFolder+'/legend_nocloud_'+exportClassifName]:
        os.remove(f)
    
    #Import the raster with the combined classifications
    baseClassif = gdal.Open(os.path.join(classifFolder,exportClassifName))
    #Mask the combined raster with the srtm data if needed and demultiply the classification based on the elevation breaks
    if not srtmMask:
        processed = new_raster_from_base(baseClassif, tempFolder+'/'+'masked0.tif', 'GTiff', -32768, gdal.GDT_Int16, bands=1)
        processed.GetRasterBand(1).WriteArray(baseClassif.GetRasterBand(1).readAsArray())
        
        #Number of pieces of the classification raster to tabulate
        pieces = 1
        
        baseClassif = None
        processed = None
    else:
        #Import the elevation raster
        elev = gdal.Open(srtmRaster)
        
        nodata = elev.GetRasterBand(1).GetNoDataValue()
        
        #Adapt the resolution of the smooth images to the baseline if needed
        elev = warp_raster(elev, baseClassif, resampleOption='average', outputURI=None, outFormat='MEM')
        
        #Number of pieces of the classification raster to tabulate
        pieces = len(breakValues) + 1
        
        #Complete the breakValues
        breakValues.insert(0, 0)
        breakValues.append(10000)
        
        #Get the size of the rasters to identify the limits of the blocks to loop through
        band = baseClassif.GetRasterBand(1)
        #Get the size of the raster
        xsize = band.XSize
        ysize = band.YSize
        #Set the block size
        BlockXSize = 256
        BlockYSize = 256
        #Get the number of blocks in x and y directions
        xBlocks = int(round(xsize/BlockXSize)) + 1
        yBlocks = int(round(ysize/BlockYSize)) + 1
        
        processed = []
        for i in range(pieces):
            #Prepare empty rasters to hold the masked rasters
            processed.append(new_raster_from_base(baseClassif, tempFolder+'/'+'masked'+str(i)+'.tif', 'GTiff', -32768, gdal.GDT_Int16, bands=1))
        
        for xStep in range(xBlocks):
            for yStep in range(yBlocks):
                
                #Import the block of the classification
                baseBlock = readRasterBlock(baseClassif, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize)
                #Import the block of the elevation raster
                elevBlock = readRasterBlock(elev, xStep*BlockXSize, yStep*BlockYSize, BlockXSize, BlockYSize)
                
                for i in range(pieces):
                    #Mask the values based on the elevation raster
                    if not nodata:
                        outArray[np.logical_or(np.isnan(elevBlock),np.logical_or(elevBlock<breakValues[i],elevBlock>=breakValues[i+1]))] = -32768
                    else:
                        outArray[np.logical_or(np.logical_or(np.isnan(elevBlock),elevBlock==nodata),np.logical_or(elevBlock<breakValues[i],elevBlock>=breakValues[i+1]))] = -32768                    
                    
                    #Write to output raster
                    processed[i].GetRasterBand(1).WriteArray(outArray, xStep*BlockXSize, yStep*BlockYSize)
        
        #Close the rasters
        for p in processed:
            processed[p] = None
        baseClassif = None
        elev = None
    '''
    #Import the grid
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(gridShape, 1) # 0 means read-only. 1 means writeable.
    # Check if shapefile is correctly loaded.
    if dataSource is None: 
        do_message("ERROR","PROCESS ABORTED: Error (1) opening shapefile " + gridShape)
        return False
    # Create layer
    gridLayer = dataSource.GetLayer(0)
    '''
    
    tabulated = []
    for i in len(pieces):
        #Tabulate area by grid cell for each of the combined classifications
        stat = img.tabulate_area(regionsFileName=gridShape, rasterFileName=os.path.join(tempFolder,'masked'+str(i)+'.tif'), 
                             bands=1, polIdFieldName=idField, numPix=False, outTxt=None, prefix='cat', numOutDecs=0, alltouch=False)
        #Add the tabulated values to the result list
        #The result is a dictionary of dictionaries. First dictionary are the polygon ids, the second [value in raster: count]
        tabulated.append(stat[0])
        
    #Get all the polygon ids for each elevation zone
    polygons = []
    #Get all the categories
    categories = set()
    for piece in tabulated:
        polygons.append(sorted(piece.keys()))
        for k in piece.keys():
            categories.update(k.keys())
    categories = list(categories)
    categories.sort()
    
    #Create empty arrays and fill them with the values
    tabArray = []
    for l in range(len(polygons)):
        tabArray.append(np.zeros((len(polygons),len(categories))))
        
        #For each polygon cell and category, get the associated count
        for p in range(len(polygons[l])):
            for c in range(len(categories)):
                try:
                    tabArray[l][p,c] = tabulated[l][polygons[l][p]][categories[c]]
                except:
                    continue
    
    #Import the confusion matrices
    
    
##########################################################################################################################
##########################################################################################################################
##FUNCTION    
def combine_classifications(classifFolder, orderClassif, legendFolder, legendPrefix, legendExt, 
                            legendDelim, tempFolder, exportClassifName, toNa=[], valuesNotInLegend=[]):
    
    #Check the computer platform
    ordi = checkPlatform()
    
    if ordi == "Windows":
        import gdal_merge as gm
    
    #Import the classifications names along with the order
    classifNames=[]
    priorities=[]
    with open(orderClassif,'r') as f:
        next(f) # skip headings
        reader=csv.reader(f,delimiter='\t')
        for priority,classif in reader:
            priorities.append(int(priority))
            classifNames.append(classif)
    
    #Remove the raster extensions if there are any
    classifNames = [re.sub('.(tif|TIF)$','',c) for c in classifNames]
    
    #Order the classifs in the correct order
    classifNames = [c for (p,c) in sorted(zip(priorities,classifNames), key=lambda pair: pair[0])]
    
    #Import all the legends into a list of dictionaries
    legends = []
    newLegends = []
    for nm in classifNames:
        dict = {} #Will hold the legend values in format Value:Classname
        with open(os.path.join(legendFolder,legendPrefix+nm+legendExt),'r') as f:
            next(f) # skip headings
            reader=csv.reader(f,delimiter=legendDelim)
            for v,c in reader:
                dict[int(v)] = re.sub('[0-9]*$','',c)
        #Add to the list
        legends.append(dict)
    
    #Establish a legend combining the info for all the classifications
    #Get the unique class names
    commonLegend = {v for l in legends for v in l.values()}
    commonLegend = {v:[] for v in commonLegend}
    #Add the values corresponding to each class
    i = 1
    for l in legends:
        for k,v in l.iteritems():
            commonLegend[v].append(k+100*i)
        i += 1 #The values of the rasters for each classification will be in a different 100s.
    
    #Loop through the classifications, and change the values to the combined legend
    i = 1
    for nm,legend in zip(classifNames,legends):
        
        if (os.path.isfile(os.path.join(classifFolder, nm+'.tif'))):
            tif = '.tif'
        else:
            tif = '.TIF'
        
        try:
            classifRaster = gdal.Open(os.path.join(classifFolder, nm+tif))
        except RuntimeError, e: 
            do_message("ERROR","PROCESS ABORTED: Error (2) opening raster file " + nm+tif + "\n" + str(e)) 
            return False
        
        if classifRaster is None: 
            do_message("ERROR","PROCESS ABORTED: Error (3) opening raster file " + nm+tif)
            return False
        
        #Get the raster band
        classifBand = classifRaster.GetRasterBand(1)
        
        #Get the no data value if any
        nodata = int(classifBand.GetNoDataValue())
        
        #Transform the band into a numpy array
        band = classifBand.ReadAsArray()
        band = band.astype(int)
        
        #Get the values in the raster in sorted order
        #(Get them all otherwise throws an error)
        legendValues = legend.keys()
        palette = list(legendValues)
        if nodata:
            palette.append(nodata)
        if valuesNotInLegend:
            for v in valuesNotInLegend:
                palette.append(v)
        palette = sorted(palette)
        
        #Create the key that gives the new values the palette should be mapped to
        mapkey = []
        for v in palette:
            if v in legendValues:
                #Get the class name for that value
                c = legend[v]
                if c in toNa:
                    mapkey.append(255)
                else:
                    mapkey.append(v+100*i)
            else:
                mapkey.append(255)
        i += 1 #The values of the rasters for each classification will be in a different 100s.
        mapkey = np.array(mapkey)
        
        #Change the values in the raster
        index = np.digitize(band, palette, right=True)
        band = mapkey[index]
        
        #Create an empty raster to export the classification that was just mapped
        outR = new_raster_from_base(classifRaster, tempFolder+'/reclass_'+nm+'.tif', 'GTiff', 255, gdal.GDT_UInt32, bands=1)
        
        #Write the re-classified values to the empty raster
        outR.GetRasterBand(1).WriteArray(band)
        
        #Close the rasters
        outR.FlushCache()
        outR = None
    
    if os.path.isfile(os.path.join(classifFolder,exportClassifName)):
        os.remove(os.path.join(classifFolder,exportClassifName))
    
    #Merge the exported rasters 
    #The processing command for merge does not allow to declare missing values...
    #processing.runalg("gdalogr:merge", {"INPUT":";".join(tarNames), "PCT":False, 
    #        "SEPARATE":True, "RTYPE":3, "OUTPUT":os.path.join(out,outName)})
    #Prepare the arguments for gdal_merge.py -- They all need to be strings
    args = ['-n', '255', '-a_nodata', '255', '-of', 'GTiff', '-ot', 'UInt32', '-co', 'COMPRESS:LZW', '-o', classifFolder+'/'+exportClassifName]
    args.extend([tempFolder+'/reclass_'+nm+'.tif' for nm in reversed(classifNames)])
    
    if ordi == "Linux":
        args = ['gdal_merge.py','-of', "GTiff"]+args
        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        p.communicate()[0] #Waits until the subprocess is finished before running the rest of the code
    elif ordi == "Windows":
        #Windows makes an error when calling gdal_merge directly, 
        #So need to specify python for the executable and to make it run in the shell
        sys.argv = args
        gm.main()
    
    #Export the combined legend to a txt file
    with open(os.path.join(legendFolder,legendPrefix+exportClassifName), "wb") as f:
        f.write('%s\t%s\n' % ('code', 'classname'))
        for k,v in commonLegend.iteritems():
            f.write('%s\t%s\n' % (v,k))
    
    #Clean the temp folder
    for nm in classifNames:
        os.remove(tempFolder+'/reclass_'+nm+'.tif')
    
    
def new_raster_from_base(base, outputURI, format, nodata, datatype, bands=None):
    """
    ---------------------------------------------------------------------------------------------
    Function : Create an empty copy of a raster from an existing one
               
    Arguments : 

    - base: gdal raster layer
        Name of the variable with the input raster to copy
    
    - outputURI: string
        Address + name of the output raster (extension should agree with format, none for memory)
        
    - format: string
        Format for the dataset (e.g. "GTiff", "MEM") 
    
    - nodata: int/float
        No data value (type should agree with raster type)
    
    - datatype: gdal data type (e.g. gdal.GDT_Int32)
        Data type for the raster
    
    - bands: [optional] int
        Number of bands for the output raster. 
        If not specified, will use the number of bands of the input raster
    
    Return value : A gdal raster variable filled with the nodata value
    ---------------------------------------------------------------------------------------------
    """    
    
    cols = base.RasterXSize
    rows = base.RasterYSize
    projection = base.GetProjection()
    geotransform = base.GetGeoTransform()
    if not bands:
        bands = base.RasterCount

    driver = gdal.GetDriverByName(format)
    
    if format == "GTiff":
        new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype, options=['COMPRESS=LZW'])
    else:
        new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype)
    new_raster.SetProjection(projection)
    new_raster.SetGeoTransform(geotransform)

    for i in range(bands):
        new_raster.GetRasterBand(i + 1).SetNoDataValue(nodata)
        new_raster.GetRasterBand(i + 1).Fill(nodata)

    return new_raster

def do_message(msgType, msgStr):
    """
    ---------------------------------------------------------------------------------------------
    Function : prints an error message
    ---------------------------------------------------------------------------------------------
    """
    print "\n============================================================"
    print msgType + ": " + msgStr
    print "============================================================"    


def checkPlatform():
    '''
    Function : Checks for the current platform on which the code is running
    '''
    
    currentPlatform = platform.system()
    
    if not currentPlatform == "Linux" and not currentPlatform == "Windows":
        return
    
    return(currentPlatform)

def readRasterBlock(src, xStart, yStart, xBlockSize, yBlockSize, band=1):
    '''
    Function to read a block of data from a gdal Dataset. It will return the block as a numpy array
    
    src (gdal Dataset): raster to extract the block from
    xStart (int): X pixel value from which to start extraction (upper left corner)
    yStart (int): Y pixel value from which to start extraction (upper left corner)
    xBlockSize (int): X size of the block to extract. If the block is larger than the number 
            of pixels from xStart, will extract up to the raster limit
    yBlockSize (int): Y size of the block to extract. If the block is larger than the number 
            of pixels from xStart, will extract up to the raster limit
    band (int): band from which to extract the block. 1 by default.
    '''
    
    #Get the wanted band
    band = src.GetRasterBand(band)
    
    #Get the size of the raster
    xsize = band.XSize
    ysize = band.YSize
    
    if yStart + yBlockSize < ysize:
        rows = yBlockSize
    else:
        rows = ysize - yStart
    
    if xStart + xBlockSize < xsize:
        cols = xBlockSize
    else:
        cols = xsize - xStart
            
    outArray = band.ReadAsArray(xStart, yStart, cols, rows)
    
    return outArray

def warp_raster(src, dst, resampleOption='nearest', outputURI=None, outFormat='MEM'):
    """
    ---------------------------------------------------------------------------------------------
    Function : Warp a source raster to the resolution, extent and projection of a destination raster.
           
            The function returns the resulting raster. If outFormat is different from 'MEM', the 
            raster is also saved to disk using the information provided in outputURI.
            
            Inputs
            --src (gdal Dataset): source raster to be warped
            --dst (gdal Dataset): destination raster that will provide the resolution, extent and projection
            --resampleOption (string): One of 'nearest', 'bilinear', 'cubic', 'cubic spline', 'lanczos', 
                    'average', or 'mode'. Method to use to resample the pixels of the source raster
            --outputURI (string, optional): Full address and name of the output raster. If outFormat is 'MEM',
                    this argument is ignored and the function simply produces a raster in memory. 
                    The extension for the output file should match the outFormat.
            
            --outFormat (string, optional): Format to use for the output raster from the function. 
                    Use 'GTiff' for a .tif output. Default creates a raster in memory.
    ---------------------------------------------------------------------------------------------
    """  
    if not type(src) is gdal.Dataset:
        return False
    
    if not type(dst) is gdal.Dataset:
        return False
    
    #Define resampling options
    resampleOptions = {'nearest': gdalconst.GRA_NearestNeighbour, 'bilinear':gdalconst.GRA_Bilinear, 
                   'cubic':gdalconst.GRA_Cubic, 'cubic spline':gdalconst.GRA_CubicSpline, 
                   'lanczos':gdalconst.GRA_Lanczos, 'average':gdalconst.GRA_Average, 
                   'mode':gdalconst.GRA_Mode} 
    
    if not resampleOption in resampleOptions.keys():
        return False
    
    #Raster to host the warped output
    if outFormat == 'MEM':
        rOut = new_raster_from_base(dst, 'temp', 'MEM', 0, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    else:
        rOut = new_raster_from_base(dst, outputURI, outFormat, 0, src.GetRasterBand(1).DataType, bands=src.RasterCount)
    
    #Warp: the parameters are source raster, destination raster, source projection, destination projection, resampling option 
    gdal.ReprojectImage(src, rOut, src.GetProjection(), rOut.GetProjection(), resampleOptions[resampleOption])
    
    return rOut