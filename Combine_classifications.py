from osgeo import ogr, gdal, osr, gdalconst
import os, re, csv, platform, sys
import numpy as np
import subprocess

def run_script(iface):
    #Address of the folder with the classifications
    classifFolder = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/Vietnam/LD/classifs'
    
    #Address of the folder with the legends
    legendFolder = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/Vietnam/LD/classifs/Legend'
    
    #Address of the folder for the temporary files
    tempFolder = '/media/olivier/olivier_ext/gedata_current/temp'
    
    #Legend file names prefix to add to the classification names
    legendPrefix = 'legend_'
    #Legend file extensions: could be '.csv' or '.txt'
    legendExt = '.csv'
    #Legend file delimiter
    legendDel = ','
    
    #Classes to pass to NA values in combining the classifications
    #Leave an empty list if None
    #These names should not include any trailing numbers if several categories have been made (e.g. Cloud1, Cloud2...)
    toNa = ['cloud','shadow']
    valuesNotInLegend = [0]
    
    #TAB delimited file with the priority order in which to combine the classifications
    #The file should have two columns, the first is a priority number, the second the file name of the classification
    orderClassif = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/Vietnam/LD/classifs/classif_priorities_LD.txt'
    
    #Ouput name and address for the combined files and legend
    exportClassifName = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/Vietnam/LD/classifs/combined_classif_LD.tif'
    exportLegendName = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/Vietnam/LD/classifs/legend_combined_classif_LD.txt'
    
    
    
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
    for nm in classifNames:
        dict = {} #Will hold the legend values in format Value:Classname
        with open(os.path.join(legendFolder,legendPrefix+nm+legendExt),'r') as f:
            next(f) # skip headings
            reader=csv.reader(f,delimiter=legendDel)
            for v,c in reader:
                dict[int(v)] = re.sub('[0-9]*$','',c.lower())
        #Add to the list
        legends.append(dict)
    
    #Establish a common legend for all the classifications
    #Get the unique class names and remove any of the categories being mapped to NA
    commonLegend = {v for l in legends for v in l.values() if v not in toNa}
    commonLegend = list(commonLegend)
    #Add a raster value by alphabetical order
    commonLegend = {nm:i for i,nm in enumerate(sorted(commonLegend),1)}
    
    #Loop through the classifications, remove the clouds and shadows, and change the values to the common legend
    
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
        key = []
        for v in palette:
            if v in legendValues:
                #Get the class name for that value
                c = legend[v]
                if c in toNa:
                    key.append(255)
                else:
                    key.append(commonLegend[c])
            else:
                key.append(255)
        key = np.array(key)
        
        #Change the values in the raster
        index = np.digitize(band, palette, right=True)
        band = key[index]
        
        #Create an empty raster to export the classification that was just mapped
        outR = new_raster_from_base(classifRaster, tempFolder+'/reclass_'+nm+'.tif', 'GTiff', 255, gdal.GDT_Byte, bands=1)
        
        #Write the re-classified values to the empty raster
        outR.GetRasterBand(1).WriteArray(band)
        
        #Close the rasters
        outR.FlushCache()
        outR = None
    
    if os.path.isfile(exportClassifName):
        os.remove(exportClassifName)
    
    #Stack the exported rasters 
    #The processing command for merge does not allow to declare missing values...
    #processing.runalg("gdalogr:merge", {"INPUT":";".join(tarNames), "PCT":False, 
    #        "SEPARATE":True, "RTYPE":3, "OUTPUT":os.path.join(out,outName)})
    #Prepare the arguments for gdal_merge.py -- They all need to be strings
    args = ['-n', '255', '-a_nodata', '255', '-of', 'GTiff', '-ot', 'Byte', '-co', 'COMPRESS:LZW', '-o', exportClassifName]
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
    
    #Export the common legend to a txt file
    with open(exportLegendName, "wb") as f:
        f.write('%s\t%s\n' % ('code', 'classname'))
        for k,v in commonLegend.iteritems():
            f.write('%s\t%s\n' % (v,k))
    
    #Clean the temp folder
    for nm in classifNames:
        os.remove(tempFolder+'/reclass_'+nm+'.tif')
    
    
def new_raster_from_base(base, outputURI, formatD, nodata, datatype, bands=None):
    """
    ---------------------------------------------------------------------------------------------
    Function : Create an empty copy of a raster from an existing one
               
    Arguments : 

    - base: gdal raster layer
        Name of the variable with the input raster to copy
    
    - outputURI: string
        Address + name of the output raster (extension should agree with format, none for memory)
        
    - formatD: string
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

    driver = gdal.GetDriverByName(formatD)
    
    if formatD == "GTiff":
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