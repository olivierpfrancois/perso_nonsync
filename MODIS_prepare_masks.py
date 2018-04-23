import gedata_tbox_raster as rt
import gedata_tbox_classif as cl
import os, csv, re
from osgeo import gdal

root = '/media/olivier/olivier_ext1/gedata_current/jde_coffee'
dataDir = 'MODIS/collection6/terra/Tea/data'
modisDir = 'MODIS/collection6/terra/Tea'
regions = ['SLK']
tempFolder = '/media/olivier/olivier_ext1/gedata_current/temp'

priorities = [os.path.join(root,dataDir,'SLK','classifs/classif_priorities_SLK.txt')]
classifFolders = [os.path.join(root,dataDir,'SLK','classifs')]
legendFolders = [os.path.join(root,dataDir,'SLK','classifs/Legend')]
combinedClassifNms = [os.path.join(root,dataDir,'SLK','classifs/combined_classif_SLK.tif')]
combinedLegendNms = [os.path.join(root,dataDir,'SLK','classifs/legend_combined_classif_SLK.txt')]
passToMiss = [['cloud','shadow']]

#Categories to mask
categToMask = [['tea']]

#master MODIS file to use for the mask
masterModis = [os.path.join(root,modisDir,'SLK','raw_data/MOD13Q1_2016-01-01_h25-26v8_250m_16_days_NDVI.tif')]

#Output names for the masks
masksOutNms = [os.path.join(root,modisDir,'SLK','masks/SLK_densities_tea_from_classifications.tif')]

for r in range(len(regions)):
    '''
    #Combine the classifications
    cl.combineClassif(
        orderClassif=priorities[r], 
        classifFolder=classifFolders[r], 
        legendFolder=legendFolders[r], 
        tempFolder=tempFolder, 
        outClassifName=combinedClassifNms[r], 
        outLegendName=combinedLegendNms[r], 
        legendPrefix='legend_', 
        legendExt='.csv',
        legendDel=',', 
        toNa=passToMiss[r])
    '''
    #Import the combined legend
    legendC = {}
    with open(combinedLegendNms[r],'r') as f:
        next(f) # skip headings
        reader=csv.reader(f,delimiter='\t')
        for nb,categ in reader:
            legendC[categ] = int(nb)
    
    #Import the comcbined classification
    combi = gdal.Open(combinedClassifNms[r])
    
    #Create the mask
    mask = rt.createMask(src=combi, 
                         values=[legendC[nm] for nm in categToMask[r]], 
                         masking='keep', 
                         band=1, 
                         outFormat='MEM', 
                         outR=None)
    #Close the raster
    combi = None
    
    #Import the model raster
    masterM = gdal.Open(masterModis[r])
    
    #Warp to the MODIS master at same resolution
    outMask = rt.warpRaster(src=mask, dst=masterM, 
                  resampleOption='average', 
                  outDataType=gdal.GDT_Float32,
                  outputURI=re.sub('\.tif$','_250m.tif',masksOutNms[r]), 
                  outFormat='GTiff')
    #Close the raster
    masterM = None
    
    #Downsample the raster to 1km resolution
    master1km = rt.resampleRaster(src=outMask, 
                                  downFactor=4, 
                                  resampleOption='nearest', 
                                  outputURI='', 
                                  formatR='MEM', 
                                  nodata=None, 
                                  datatype=None)
    #Close the raster
    outMask = None
    
    
    #Warp to the low res master at same resolution
    outMask = rt.warpRaster(src=mask, dst=master1km, 
                  resampleOption='average', 
                  outDataType=gdal.GDT_Float32,
                  outputURI=re.sub('\.tif$','_1km.tif',masksOutNms[r]), 
                  outFormat='GTiff')
    
    #Close the rasters
    mask = None
    master1km = None
    outMask = None