from osgeo import gdal
import os, re, csv, platform, sys
import numpy as np
import subprocess
import gedata_tbox_raster as rt

def combineClassif(orderClassif, classifFolder, legendFolder, tempFolder, 
                   outClassifName, outLegendName, 
                   legendPrefix='legend_', legendExt='.txt', legendDel='\t', 
                   toNa=[]):
    '''
    Combine the classifications listed in the orderClassif file,
    with the priorities listed in the file, and using the legends 
    in the legend folder provided.
    The legends should have the same name as the classifications, 
    save for a potential prefix to be provided.
    The output is a unified classification, with a unified legend, 
    where classes with the same (exact) name have been combined.
    All values that are not in the legends will be passed to missing.
    
    orderClassif (str): Address of the file with the classifications 
        names and the priorities for combining them. 
        The file should be .txt and tab delimited.
        It should have two columns, the first is the priority number and
        the second is the base name of the classif, without extension.
        The first row is considered as the column names.
    classifFolder (str): Address of the folder containing the 
        classifications to combine.
    legendFolder (str): Address of the folder containing the
        legends of the classifications. The legend files should
        have the same name as the classifications, save for
        a potential prefix.
        CAPS do not matter, all category names will be passed 
        to lower case.
    tempFolder (str): Address of the temporary folder where
        intermediate classifications can be saved during processing.
    outClassifName (str): Full address of the output combined 
        classification.
        The function uses the value 255 as the no data value for the 
        output.
    outLegendName (str): Full address of the output legend. 
        Should be a .txt.
    legendPrefix (str): prefix used to differentiate the legends 
        names from the classifications names.
    legendExt (str): Extension of the legend files.
    legendDel (str): Delimiter of the legend files.
    toNa (list of str): Categories in the classifications 
        to pass to no data in the output classification.
        The names provided should match exactly the information
        in the legend files.
        These names should not include any trailing numbers 
        if several categories have been made (e.g. Cloud for Cloud1, Cloud2...)
    '''
    
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
        dic = {} #Will hold the legend values in format Value:Classname
        with open(os.path.join(legendFolder,legendPrefix+nm+legendExt),'r') as f:
            next(f) # skip headings
            reader=csv.reader(f,delimiter=legendDel)
            for v,c in reader:
                dic[int(v)] = re.sub('[0-9]*$','',c.lower())
        #Add to the list
        legends.append(dic)
    
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
        
        #Transform the band into a numpy array
        band = classifBand.ReadAsArray()
        band = band.astype(int)
        
        #Map the old values to the new common values
        for v in np.unique(band):
            if v in legend:
                #Get the class name for that value
                c = legend[v]
                if c in toNa:
                    replaceVal = 255
                else:
                    replaceVal = commonLegend[c]
            else:
                replaceVal = 255
            
            band[band==v] = replaceVal
        
        #Create an empty raster to export the classification that was just mapped
        outR = rt.newRasterFromBase(classifRaster, tempFolder+'/reclass_'+nm+'.tif', 'GTiff', 255, gdal.GDT_Byte, bands=1)
        
        #Write the re-classified values to the empty raster
        outR.GetRasterBand(1).WriteArray(band)
        
        #Close the rasters
        outR.FlushCache()
        outR = None
    
    if os.path.isfile(outClassifName):
        os.remove(outClassifName)
    
    #Stack the exported rasters 
    #Prepare the arguments for gdal_merge.py -- They all need to be strings
    args = ['-n', '255', '-a_nodata', '255', '-of', 'GTiff', '-ot', 'Byte', '-co', 'COMPRESS:LZW', '-o', outClassifName]
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
    #Get the categories in order
    outLegend = {}
    outCat = []
    for k,v in commonLegend.iteritems():
        outLegend[v] = k
        outCat.append(v)
    outCat.sort()
    with open(outLegendName, "wb") as f:
        f.write('%s\t%s\n' % ('code', 'classname'))
        for c in outCat:
            f.write('%s\t%s\n' % (c,outLegend[c]))
    
    #Clean the temp folder
    for nm in classifNames:
        os.remove(tempFolder+'/reclass_'+nm+'.tif')

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


if __name__=='__main__':
    
    root = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6/terra/Tea/data/SLK/classifs'
    
    combineClassif(
        orderClassif=root+'/classif_priorities_SLK.txt', 
        classifFolder=root, 
        legendFolder=root+'/Legend', 
        tempFolder='/media/olivier/olivier_ext1/gedata_current/temp', 
        outClassifName=root+'/combined_classif_SLK.tif', 
        outLegendName=root+'/legend_combined_classif_SLK.txt', 
        legendPrefix='legend_', 
        legendExt='.csv',
        legendDel=',', 
        toNa=['cloud','shadow'])
    
    