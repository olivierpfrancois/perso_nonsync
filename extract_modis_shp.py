# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" Your Description of the script goes here """

import sys
# Block python from writing pyc files
sys.dont_write_bytecode = True

# Imports
from datetime import datetime, timedelta
from osgeo import gdal, ogr
import os, re, osr
import numpy as np
from csv import DictWriter
import gedata_tbox_raster as rt

sys.path.append('/home/olivier/OlivierGithub/QGIS-scripts') #('C:/Users/Olivier/OlivierGithub')
import projection_function as proj
import image_proc_functions as img


def main():
    
    #Root directory for the modis data
    root = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6/Brazil' #'E:/gedata_current/jde_coffee/MODIS/collection6'
    
    #Full pqth to folder to use for temporary files (should be empty)
    tempDir = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/Temp'
    
    #Destination directory of the extracted modis information
    dst = root
    
    ##REGIONS parameters
    #Regions to process inputs    !!!!SHOULD BE IN EPSG 4326 PROJECTION
    #Names of the regions (also folders names) 
    states = ["CER","CHA","CO","ES","MO","SDM","SP","ZM"]
    #Varieties in each case
    varieties = [['arabica'],['arabica'],['arabica'],['arabica','robusta'],['arabica'],['arabica'],['arabica'],['arabica','robusta']]
    #Name of subfolder where the smoothed images are located
    statesSmoothFolder = 'smooth_data'
    
    #Address of the shapefile with the polygons in which to extract
    shp = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/data/Brazil/brazil_admin_boundaries/microregions.shp' # br_municipalities.shp
    #Name of the field with the unique ids of the polygons
    idAttribute = 'MICRO_REGI' #'CD_GEOCMU'
    
    #Start date for extraction
    startd = '2006-01-01'
    #End date for extraction
    endd = '2018-01-01'
    
    masks = [['masks/CER_densities_arabica_from_classifications_250m.tif'],
             ['masks/CHA_densities_arabica_from_classifications_250m.tif'],
             ['masks/CO_densities_arabica_from_classifications_250m.tif'],
             ['masks/ES_densities_arabica_from_classifications_250m.tif','masks/ES_densities_robusta_from_classifications_250m.tif'],
             ['masks/MO_densities_arabica_from_classifications_250m.tif'],
             ['masks/SDM_densities_arabica_from_classifications_250m.tif'],
             ['masks/SP_densities_arabica_from_classifications_250m.tif'],
             ['masks/ZM_densities_arabica_from_classifications_250m.tif','masks/ZM_densities_robusta_from_classifications_250m.tif']]
    
    exportFormat = 'long'
    
    extractModis(root=root, regions=states, varieties=varieties, regionsIn=statesSmoothFolder, outFolder=dst, masks=masks, 
                 startExtract=startd, endExtract=endd, shpName=shp, attr=idAttribute, tempDir=tempDir, exportFormat=exportFormat, epsgShp=4674)



def extractModis(root, regions, varieties, regionsIn, 
                 outFolder, masks, startExtract, endExtract, 
                 shpName, attr, tempDir, exportFormat, 
                 epsgShp=None):
    
    #Transform into date format
    if not endExtract:
        endExtract = datetime.now()
        endExtract = endExtract.date()
    else:
        endExtract = datetime.strptime(endExtract, '%Y-%m-%d').date()
    if not startExtract:
        startExtract = endExtract - timedelta(days=30)
    else:
        startExtract = datetime.strptime(startExtract, '%Y-%m-%d').date()
    
    
    #Open the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shpName, 0) # 0 means read-only. 1 means writeable.

    # Check if shapefile is correctly loaded.
    if dataSource is None: 
        print("ERROR","PROCESS ABORTED: Error (1) opening shapefile " + shpName)
        return False
    
    # Get layer
    inVecLayer = dataSource.GetLayer(0)
    
    #Get projection if possible
    inEpsg = proj.getProjection(dataSource)
    if not inEpsg and not epsgShp:
        print('Could not identify projection for the shapefile. Specify a value for epsgShp parameter')
        return False
    elif not inEpsg:
        inEpsg = epsgShp
    
    if not isinstance(inEpsg,int):
        inEpsg = int(inEpsg)
    
    if not inEpsg == 4326:
        
        #Create the transform
        prj1 = osr.SpatialReference()
        prj1.ImportFromEPSG(inEpsg)
        prj2 = osr.SpatialReference()
        prj2.ImportFromEPSG(4326)
        coordTrans = osr.CoordinateTransformation(prj1,prj2)
        
        #create output datasource
        outdriver = ogr.GetDriverByName('ESRI Shapefile')
        srcOut = outdriver.CreateDataSource(os.path.join(tempDir,'temp.shp'))
        #Create a new layer
        lyrOut = srcOut.CreateLayer("layer 1", srs=prj2, geom_type=inVecLayer.GetGeomType())
        
        #Add fields
        lyr1Defn = inVecLayer.GetLayerDefn()
        for i in range(0, lyr1Defn.GetFieldCount()):
            fieldDefn = lyr1Defn.GetFieldDefn(i)
            lyrOut.CreateField(fieldDefn)
        
        # get the output layer's feature definition
        lyrOutDefn = lyrOut.GetLayerDefn()
        
        # loop through the input features
        inFeature = inVecLayer.GetNextFeature()
        while inFeature:
            # get the input geometry
            geom = inFeature.GetGeometryRef()
            # reproject the geometry
            geom.Transform(coordTrans)
            # create a new feature
            outFeature = ogr.Feature(lyrOutDefn)
            # set the geometry and attribute
            outFeature.SetGeometry(geom)
            for i in range(0, lyrOutDefn.GetFieldCount()):
                outFeature.SetField(lyrOutDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
            # add the feature to the shapefile
            lyrOut.CreateFeature(outFeature)
            # dereference the features and get the next input feature
            outFeature = None
            inFeature = inVecLayer.GetNextFeature()
        
        #Close the two shapefiles
        dataSource.Destroy()
        inVecLayer = None
        srcOut.Destroy()
        lyrOut = None
    
    #Empty dictionary to hold the results.
    #Every municipality covered by one of the regions will have a coffee weight
    coffeeWeights = {'arabica':{}, 'robusta':{}}
    #Every municipality will have a total number of modis pixels
    totPixels = {}
    #Every municipality covered by one of the regions will have a share of pixels with coffee in it
    coffeeShare = {'arabica':{}, 'robusta':{}}
    #Every municipality covered by one of the regions will have a coffee count
    #Every date will have a dictionary containing the municipalities and values
    imgStats = {'arabica':{}, 'robusta':{}}
    
    
    #Loop through the regions to extract the modis data for each
    for r in range(len(regions)):
        print('Extracting region '+str(regions[r])+'...')
        
        #Import all the smooth modis images on disk
        onDisk = [f for f in os.listdir(os.path.join(root,regions[r],regionsIn)) if f.endswith('.tif')]
        
        #Dates of these files
        datesAll = [re.search('_([0-9]{4}-[0-9]{2}-[0-9]{2})', f).group(1) for f in onDisk]
        #Transform into date format
        datesAll = [datetime.strptime(d, '%Y-%m-%d').date() for d in datesAll]
            
        #Keep only the files and dates within the dates to process
        onDisk = [f for f,d in zip(onDisk,datesAll) if d >= startExtract and d <= endExtract]
        datesAll = [d for d in datesAll if d >= startExtract and d <= endExtract]
        #Sort the two list by date
        datesAll, onDisk = (list(x) for x in zip(*sorted(zip(datesAll, onDisk))))
        
        if not onDisk:
            print('no modis images to process in '+regions[r])
            continue
        
        #Get a base image as a reference for format and extent
        baseImg = gdal.Open(os.path.join(root,regions[r],regionsIn,onDisk[0]))
        
        #Get the corner latitude and longitude for the raster
        cols = baseImg.RasterXSize # width
        rows = baseImg.RasterYSize # height
        gt = baseImg.GetGeoTransform()
        minLonR = gt[0]
        minLatR = gt[3] + cols*gt[4] + rows*gt[5] 
        maxLonR = gt[0] + cols*gt[1] + rows*gt[2]
        maxLatR = gt[3]
        
        #Close the raster
        baseImg = None
        
        #Crop the shapefile by the extent of the raster
        inShp = tempDir+'/temp.shp'
        outShp = tempDir+'/temp_crop.shp'
        #Call ogr2ogr using the os system command
        cmd = 'ogr2ogr -f "ESRI Shapefile" %s %s -clipsrc %f %f %f %f -overwrite' % (outShp, inShp, minLonR, minLatR, maxLonR, maxLatR)
        os.system(cmd)
        
        
        ###Extract the number of modis pixels per cropped municipalities
        #Import the mask
        weightsRaster = gdal.Open(root+'/'+regions[r]+'/'+masks[r][0])
        #Create an empty copy for storing the temporary weighted image
        weighted = rt.newRasterFromBase(weightsRaster, tempDir+'/temp.tif', 'GTiff', 0., gdal.GDT_Float32, bands=1)
        #Import coffee weights band as array
        weightsBand = weightsRaster.GetRasterBand(1)
        weightsArray = weightsBand.ReadAsArray()
        #Recast the type to be sure
        weightsArray = weightsArray.astype(np.float32)
        #Replace all data by 1
        weightsArray[:] = 1.
        #Export to temporary raster
        weighted.GetRasterBand(1).WriteArray(weightsArray)
        #Close the temporary raster
        weighted.FlushCache()
        weighted = None
        #Compute the count of pixels for each polygon of the cropped municipalities
        stat = img.raster_stats(regionsFileName=os.path.join(tempDir,'temp_crop.shp'), rasterFileName=tempDir+'/temp.tif',  
                                         polIdFieldName=attr, statistics=['sum'], outTxt=None, addShp=False, addSuffix='', 
                                         numOutDecs=2, alltouch=True)
        if not stat:
            print('error extracting coffee sums per polygon')
            continue
        stat = stat[0]
        #Add the counts to the total pixel counts
        for polyid in stat.keys():
            if not str(polyid) in totPixels:
                totPixels[str(polyid)] = stat[polyid]['sum']
            else:
                totPixels[str(polyid)] = totPixels[str(polyid)] + stat[polyid]['sum']
        #Remove existing temporary raster if any
        
        try:
            os.remove(os.path.join(tempDir,'temp.tif'))
        except OSError:
            pass
        
        
        
        for v in range(len(varieties[r])):
            
            for date in datesAll:
                if not date.strftime('%Y-%m-%d') in imgStats[varieties[r][v]]:
                    imgStats[varieties[r][v]][date.strftime('%Y-%m-%d')] = {}
                    imgStats[varieties[r][v]][date.strftime('%Y-%m-%d')]['date'] = date.strftime('%Y-%m-%d')
            
            ###Extract the number of pixels with some coffee per cropped municipalities
            #Import the mask
            weightsRaster = gdal.Open(root+'/'+regions[r]+'/'+masks[r][v])
            #Create an empty copy for storing the temporary weighted image
            weighted = rt.newRasterFromBase(weightsRaster, tempDir+'/temp.tif', 'GTiff', 0., gdal.GDT_Float32, bands=1)
            #Import coffee weights band as array
            weightsBand = weightsRaster.GetRasterBand(1)
            weightsArray = weightsBand.ReadAsArray()
            #Recast the type to be sure
            weightsArray = weightsArray.astype(np.float32)
            #Replace the no data value by 0 and all non zero by 1
            nodataW = weightsBand.GetNoDataValue()
            if not nodataW:
                weightsArray[np.logical_or(np.isnan(weightsArray),np.logical_or(weightsArray>1.,weightsArray<0.))] = 0
            else:
                weightsArray[np.logical_or(np.logical_or(np.isnan(weightsArray),weightsArray==nodataW),np.logical_or(weightsArray>1.,weightsArray<0.))] = 0
            nodataW = None
            #Replace all positive coffee values by 1
            weightsArray[weightsArray > 0.] = 1
            #Export to temporary raster
            weighted.GetRasterBand(1).WriteArray(weightsArray)
            #Close the temporary raster
            weighted.FlushCache()
            weighted = None
            
            '''
            if regions[r] == 'ES' and varieties[r][v] == 'robusta':
                try:
                    os.remove(os.path.join(tempDir,'tempES.tif'))
                except OSError:
                    pass
                weighted = rt.newRasterFromBase(weightsRaster, tempDir+'/tempES.tif', 'GTiff', 0., gdal.GDT_Float32, bands=1)
                #Export to temporary raster
                weighted.GetRasterBand(1).WriteArray(weightsArray)
                #Close the temporary raster
                weighted.FlushCache()
                weighted = None
            ''' 
            
            #Compute the count of pixels for each polygon of the cropped municipalities
            stat = img.raster_stats(regionsFileName=os.path.join(tempDir,'temp_crop.shp'), rasterFileName=tempDir+'/temp.tif',  
                                             polIdFieldName=attr, statistics=['sum'], outTxt=None, addShp=False, addSuffix='', 
                                             numOutDecs=2, alltouch=True)
            if not stat:
                print('error extracting coffee sums per polygon')
                continue
            stat = stat[0]
            #Add the sums to the share of coffee
            for polyid in stat.keys():
                #Do not consider the polygon if all values are missing or 0
                if np.isnan(stat[polyid]['sum']) or stat[polyid]['sum'] == 0:
                    continue
                if not str(polyid) in coffeeWeights[varieties[r][v]]:
                    coffeeShare[varieties[r][v]][str(polyid)] = stat[polyid]['sum']
                else:
                    coffeeShare[varieties[r][v]][str(polyid)] = coffeeShare[varieties[r][v]][str(polyid)] + stat[polyid]['sum']
            #Remove existing temporary raster if any
            try:
                os.remove(os.path.join(tempDir,'temp.tif'))
            except OSError:
                pass
            weightsRaster = None
            
            ###Extract the sum of coffee weights per cropped municipality
            #Import the mask
            weightsRaster = gdal.Open(root+'/'+regions[r]+'/'+masks[r][v])
            #Create an empty copy for storing the temporary weighted image
            weighted = rt.newRasterFromBase(weightsRaster, tempDir+'/temp.tif', 'GTiff', 0, gdal.GDT_Float32, bands=1)
            #Import coffee weights band as array
            weightsBand = weightsRaster.GetRasterBand(1)
            weightsArray = weightsBand.ReadAsArray()
            #Recast the type to be sure
            weightsArray = weightsArray.astype(np.float32)
            #Replace the no data values by 0
            nodataW = weightsBand.GetNoDataValue()
            if not nodataW:
                weightsArray[np.logical_or(np.isnan(weightsArray),np.logical_or(weightsArray>1.,weightsArray<0.))] = 0.
            else:
                weightsArray[np.logical_or(np.logical_or(np.isnan(weightsArray),weightsArray==nodataW),np.logical_or(weightsArray>1.,weightsArray<0.))] = 0.
            nodataW = None
            #Export to temporary raster
            weighted.GetRasterBand(1).WriteArray(weightsArray)
            #Close the temporary raster
            weighted.FlushCache()
            weighted = None
            #Compute the sum for each polygon of the coffee
            stat = img.raster_stats(regionsFileName=os.path.join(tempDir,'temp_crop.shp'), rasterFileName=tempDir+'/temp.tif',  
                                             polIdFieldName=attr, statistics=['sum'], outTxt=None, addShp=False, addSuffix='', 
                                             numOutDecs=2, alltouch=True)
            if not stat:
                print('error extracting coffee sums per polygon')
                continue
            stat = stat[0]
            #Add the weights to the coffee weights dictionary
            for polyid in stat.keys():
                #Do not consider the polygon if all values are missing or 0
                if np.isnan(stat[polyid]['sum']) or stat[polyid]['sum'] == 0:
                    continue
                if not str(polyid) in coffeeWeights[varieties[r][v]]:
                    coffeeWeights[varieties[r][v]][str(polyid)] = stat[polyid]['sum']
                else:
                    coffeeWeights[varieties[r][v]][str(polyid)] = coffeeWeights[varieties[r][v]][str(polyid)] + stat[polyid]['sum']
            
            #Loop through the images
            for i, date in zip(onDisk,datesAll):
                print(date)
                
                #Import the image
                baseImg = gdal.Open(os.path.join(root,regions[r],regionsIn,i))
                
                #Remove existing temporary raster if any
                try:
                    os.remove(os.path.join(tempDir,'temp.tif'))
                except OSError:
                    pass
                
                #Create an empty copy for storing the temporary weighted image
                weighted = rt.newRasterFromBase(baseImg, tempDir+'/temp.tif', 'GTiff', np.nan, gdal.GDT_Float32, bands=1)
                
                #Import image band as array
                imgBand = baseImg.GetRasterBand(1)
                imgArray = imgBand.ReadAsArray()
                
                #Recast the type to be sure
                imgArray = imgArray.astype(np.float32)
                
                #Replace the no data values by 0
                nodataI = imgBand.GetNoDataValue()
                if not nodataI:
                    imgArray[np.logical_or(np.isnan(imgArray),np.logical_or(imgArray>1.,imgArray<-1.))] = 0.
                else:
                    imgArray[np.logical_or(np.logical_or(np.isnan(imgArray),imgArray==nodataI),np.logical_or(imgArray>1.,imgArray<-1.))] = 0.
                
                #Change the resolution of the raster to match the images if needed
                geoSmooth = weightsRaster.GetGeoTransform()
                geoBase = baseImg.GetGeoTransform()
                reproject = [1 for a,b in zip(geoSmooth,geoBase) if not a==b]
                if reproject:
                    weightsReproj = rt.warpRaster(weightsRaster, baseImg, resampleOption='nearest', outputURI=None, outFormat='MEM')
                
                else:
                    weightsReproj = weightsRaster
                    
                #Import coffee weights band as array
                weightsBand = weightsReproj.GetRasterBand(1)
                weightsArray = weightsBand.ReadAsArray()
                
                #Recast the type to be sure
                weightsArray = weightsArray.astype(np.float32)
                
                #Replace the no data values by 0
                nodataW = weightsBand.GetNoDataValue()
                if not nodataW:
                    weightsArray[np.logical_or(np.isnan(weightsArray),np.logical_or(weightsArray>1.,weightsArray<0.))] = 0.
                else:
                    weightsArray[np.logical_or(np.logical_or(np.isnan(weightsArray),weightsArray==nodataW),np.logical_or(weightsArray>1.,weightsArray<0.))] = 0.
                
                #Estimate the weighted sum for each pixel
                weightedArray = np.multiply(weightsArray, imgArray)
                
                #Export to temporary raster
                weighted.GetRasterBand(1).WriteArray(weightedArray)
                
                #Close the temporary raster
                weighted.FlushCache()
                weighted = None
                weightsReproj = None
                
                stat = img.raster_stats(regionsFileName=os.path.join(tempDir,'temp_crop.shp'), rasterFileName=tempDir+'/temp.tif',  
                                        polIdFieldName=attr, statistics=['sum'], outTxt=None, addShp=False, addSuffix='', 
                                        numOutDecs=2, alltouch=True)
                
                if not stat:
                    print('error extracting ndvi sums per polygon for date '+date.strftime('%Y-%m-%d'))
                    continue
                
                stat = stat[0]
                #Output is a list (here one element) where each element is a dictionary. Keys are the shp ids. Each then is 
                #also a dictionary, with each key being a summary statistics
                
                for polyid in stat.keys():
                    #Do not consider the polygon if not contained in the raster bounding box or if there is no coffee
                    #if not polyid in insideFeat or np.isnan(stat[polyid]['sum']) or not str(polyid) in coffeeWeights[varieties[r][v]].keys():
                    if np.isnan(stat[polyid]['sum']) or not str(polyid) in coffeeWeights[varieties[r][v]].keys():
                        continue
                    
                    if not str(polyid) in imgStats[varieties[r][v]][date.strftime('%Y-%m-%d')]:
                        imgStats[varieties[r][v]][date.strftime('%Y-%m-%d')][str(polyid)] = stat[polyid]['sum']
                    else:
                        imgStats[varieties[r][v]][date.strftime('%Y-%m-%d')][str(polyid)] = (imgStats[varieties[r][v]][date.strftime('%Y-%m-%d')][str(polyid)] + 
                                                                                             stat[polyid]['sum'])
            
            #Remove the temporary file
            try:
                os.remove(os.path.join(tempDir,'temp.tif'))
            except OSError:
                pass
            
            #Close the mask 
            weightsRaster = None
        
        #Remove the temporary shapefile
        try:
            driver.DeleteDataSource(os.path.join(tempDir,'temp_crop.shp'))
        except OSError:
            pass
        
    #Combine the data from the coffee weights and the sum of the ndvi
    for v in ['arabica','robusta']:
        if imgStats[v]:
            for date in imgStats[v].keys():
                for polyid in imgStats[v][date].keys():
                    
                    if not polyid == 'date':
                        # Remove the region/municipality if the share of pixels with some coffee is 
                        # less than 5% of the total pixels covered by MODIS
                        if not polyid in coffeeShare[v].keys() or coffeeShare[v][polyid]/totPixels[polyid] < 0.05:
                            del imgStats[v][date][polyid]
                            continue
                        
                        imgStats[v][date][polyid] = imgStats[v][date][polyid]/coffeeWeights[v][polyid]
            
            #Export all the dates and polygon ids to a text file
            outNm = 'Weighted_avgndvi_permicroregion_'+v+'_'+startExtract.strftime('%Y-%m-%d')+'_'+endExtract.strftime('%Y-%m-%d')+'.txt'
            
            if exportFormat == 'wide':
                #order the output by polygon id in a list
                #Each element is a dictionary with the fields to be exported (polygon id, date and value)
                transposed = {}
                polyids = set()
                for date in imgStats[v].keys():
                    polyids.update(imgStats[v][date].keys())
                polyids.remove('date')
                for polyid in polyids:
                    transposed[polyid] = {}
                    transposed[polyid]['id'] = polyid
                for date in imgStats[v].keys():
                    keys = imgStats[v][date].keys()
                    keys.remove('date')
                    for polyid in keys:
                        transposed[polyid][date] = imgStats[v][date][polyid]
                imgStats[v] = transposed
                
                #order the output by date in a list. Each element is an element of the original dictionary and will be exported
                out = []
                colnames = set()
                
                for k in imgStats[v].keys():
                    out.append(imgStats[v][k])
                    colnames.update(imgStats[v][k].keys())
                colnames = list(colnames)
                if 'date' in colnames:
                    colnames.remove('date')
                    colnames.sort(key=float)
                    colnames.insert(0, 'date')
                else:
                    colnames.remove('id')
                    colnames.sort(key=lambda d: datetime.strptime(d, '%Y-%m-%d'))
                    colnames.insert(0, 'id')
            
            elif exportFormat == 'long':
                out = []
                for date in imgStats[v].keys():
                    polyids = imgStats[v][date].keys()
                    polyids.remove('date')
                    for polyid in polyids:
                        out.append({'id':polyid, 'date':date, 'value':imgStats[v][date][polyid]})
                
                colnames = ['id','date','value']
            
            with open(os.path.join(outFolder,outNm), "w") as f:
                dict_writer = DictWriter(f, colnames, extrasaction='ignore', delimiter="\t", restval="0")
                dict_writer.writeheader()
                for p in out:
                    dict_writer.writerow(p)
        
    #Remove the temporary shapefile
    try:
        driver.DeleteDataSource(os.path.join(tempDir,'temp.shp'))
    except OSError:
        pass

if __name__ == '__main__':
    main()