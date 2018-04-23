# Implements the gapfill algorithm as in the gapfill package in R

#from mpl_toolkits.basemap import Basemap
import mpl_toolkits
mpl_toolkits.__path__.append('/usr/lib/python2.7/dist-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap, pyproj
import gdal, math
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as col
import numpy as np

def mapModisRanking(mapSize, mapTitle, mapFile,
                    boundaryFile, notePosition, 
                    legendPosition, outName, outRes=200,
                    backgroundLabel='',
                    citiesFile=None, citiesField=None,
                    citiesLabelSize=None, 
                    citiesMarkerSize=None):
    
    '''
    Maps the ranking raster and export to file
    
    mapSize (tuple of num): Size of the ap in inches
    mapTitle (str): Title of the map 
    mapFile (str): Full address of the raster with 
        the ranking information
    boundaryFile (str): Full address of the shapefile 
        with the boundaries information
    notePosition (tuple of num): Two values between 0 and 
        1 giving the position of the top left corner of the 
        note on the chart. 
        The value are relative to the chart, so (0,0) is the
        bottom left corner and (1,1) is the top right corner. 
    legendPosition (tuple of num): Two values giving the 
        position of the legend. Same as notePosition
    outName (str): Full address of the output name for the
        chart (.png)
    outRes (int): Dpi resolution for the output chart
    backgroundLabel (str): Label for the white background
        pixels
    citiesFile (str): Full address of the shapefile with the 
        information on the cities.
    citiesField (str): Name of the field in the shapefile with 
        name of the cities to use as label
    citiesLabelSize (int): Size of the labels
    '''

    #Create new figure window
    fig = plt.figure(figsize=mapSize)  # a new figure window
    ax = fig.add_subplot(1, 1, 1)  # specify (nrows, ncols, axnum)
    
    #Remove frame of subplot
    ax.axis('off')
    
    #Add title
    ax.set_title(mapTitle, 
                 fontsize=24, weight='bold', y=1.02)
    
    # Read the data and metadata
    datafile = gdal.Open(mapFile)
    bnd1 = datafile.GetRasterBand(1).ReadAsArray()
    
    #Get no data value
    nodata = datafile.GetRasterBand(1).GetNoDataValue()
    
    #Change data type and remove no data value from raster
    bnd1 = bnd1.astype(float)
    bnd1[bnd1==nodata] = np.nan
    
    #Get raster size, projection, and resolution
    nx = datafile.RasterXSize # Raster xsize
    ny = datafile.RasterYSize # Raster ysize
    gt = datafile.GetGeoTransform()
    proj = datafile.GetProjection()
    xres = gt[1]
    yres = gt[5]
    
    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * nx) - xres * 0.5
    ymin = gt[3] + (yres * ny) + yres * 0.5
    ymax = gt[3] - yres * 0.5
    
    # create a grid of lat/lon coordinates in the original projection
    (lon_source,lat_source) = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
    
    #Create the basemap
    mapR = Basemap(projection='cyl',llcrnrlat=ymin,urcrnrlat=ymax,\
                llcrnrlon=xmin,urcrnrlon=xmax , resolution='i', ax=ax)
    
    #Prepare the color map
    #Colors and values for scale
    colorsScale = [(255,254,141),(239,48,166),(197,38,182),(113,28,198),(19,0,236),(28,67,198),(5,251,255)]
    stopsScale = [-3,10,30,50,70,90,110]
    #Normalize the colors to 0-1
    norm = mpl.colors.Normalize(vmin=0.,vmax=255.)
    colorsScale = [tuple(norm(v) for v in T) for T in colorsScale]
    #Normalize the values for the scale
    norm = mpl.colors.Normalize(vmin=-3.,vmax=110.)
    stopsScale = [norm(v) for v in stopsScale]
    #Combine into color scale for mapping
    purpleBlue = zip(stopsScale,colorsScale)
    #Create the segmented color map
    purpleBlueLinear = col.LinearSegmentedColormap.from_list('purpleBlue',purpleBlue, N=256, gamma=1.0)
    #Create labels for colormap
    purpleBlueLabels = ['Below min', 
                        'Min of the 10 ref. years',
                        '',
                        'Median of the 10 ref. years',
                        '',
                        'Max of the 10 ref. Years',
                        'Above max']
    
    # project in the original Basemap and plot with pcolormesh
    mapR.pcolormesh(lon_source,lat_source,bnd1.T, cmap=purpleBlueLinear)
    
    #Add boundary
    aoi_info = mapR.readshapefile(boundaryFile, 
                                  'aoi', color='black',
                                  linewidth=1.3)
    
    if citiesFile:
        #Add major cities
        cities_info = mapR.readshapefile(citiesFile, 
                                         'cities')
        
        #Add the city names as labels
        cityFont = {'fontname':'Arial', 'size':str(citiesLabelSize), 
                    'color':'black', 'weight':'bold'}
        for info, city in zip(mapR.cities_info, mapR.cities):
            mapR.plot(city[0], city[1], marker='o', 
                      color='black', markersize=citiesMarkerSize, 
                      markeredgewidth=2) #'o' for circle, '.' for point
            plt.text(city[0]+0.01, city[1]+0.005, 
                     unicode(info[citiesField], 'utf-8'), **cityFont)
    
    #Add scale bar
    scaleParam = {'startLon': xmin+(xmax-xmin)*0.62, 
                  'startLat': ymin+0.01,
                  'lengthKm': 100, 
                  'yoffset': 0.02}
    #Get initial lon lat in map units
    lon1,lat1 = mapR(scaleParam['startLon'],scaleParam['startLat'],inverse=True)
    #Get final lon lat from distance
    gc = pyproj.Geod(a=mapR.rmajor,b=mapR.rminor)
    lon2, lat2, az = gc.fwd(lon1,lat1,90,scaleParam['lengthKm']*1000)
    #Get back the final lon lat in map units
    x2,y2 = mapR(lon2,lat2,inverse=False)
    #Plot the lines for the scale
    barHeight = abs(scaleParam['startLon']-x2)/100.
    mapR.plot([scaleParam['startLon'],x2],
              [scaleParam['startLat'],scaleParam['startLat']],color='k')
    mapR.plot([scaleParam['startLon'],scaleParam['startLon']],
              [scaleParam['startLat']-barHeight,
               scaleParam['startLat']+barHeight],color='k')
    mapR.plot([x2,x2],
              [scaleParam['startLat']-barHeight,
               scaleParam['startLat']+barHeight],color='k')
    scaleFont = {'fontname':'Arial', 'size':'18', 'color':'black', 'weight':'bold', 
                 'horizontalalignment':'center'}
    plt.text(scaleParam['startLon'],scaleParam['startLat']+barHeight+0.01,
             '0', **scaleFont)
    plt.text(x2,scaleParam['startLat']+barHeight+0.01,
             '%s km' %(scaleParam['lengthKm']),    
             **scaleFont)
    
    
    #Add note
    commentFont = {'fontname':'Arial', 'size':'18', 'color':'black'}
    comment = 'Notes:'\
        '\n--Most pixels contain other land uses \nbesides coffee.'\
        '\n--The index shows health of vegetation in \neach pixel compared to reference years, \n'\
        'not coffee production directly. '\
        '\n--Vegetative health is affected mostly by \n'\
        'natural factors (rain, etc.) but can also be \naffected by human intervention (pruning, \netc.).'
    plt.text(notePosition[0], notePosition[1], comment, 
             horizontalalignment='left',
             verticalalignment='center',
             transform = ax.transAxes, 
             **commentFont)
    
    #Prepare color legend
    if backgroundLabel:
        cmap = mpl.colors.ListedColormap([(1,1,1)]+colorsScale)
        #bounds = range(len(colorsScale)+1)
    else:
        cmap = mpl.colors.ListedColormap(colorsScale)
    bounds = range(cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    #Add axis for the color legend
    legendFont = {'fontsize': 18,
                  'fontweight': 'bold',
                  'verticalalignment': 'center'}
    axColors = fig.add_axes([legendPosition[0],legendPosition[1], 0.05, 0.15]) #left, bottom, width, height] in fractions of figure width and height
    cb = mpl.colorbar.ColorbarBase(axColors, cmap=cmap,
                                    norm=norm,
                                    boundaries=bounds,
                                    ticks=[y+0.5 for y in range(len(colorsScale)+1)],
                                    spacing='uniform',
                                    orientation='vertical')
    cb.ax.set_title('Legend', fontsize=20, weight='bold', x=0.7, y=1.05)
    axColors.tick_params(axis=u'both', which=u'both',length=0)
    if backgroundLabel:
        legendLabels = [backgroundLabel]+purpleBlueLabels
    else:
        legendLabels = purpleBlueLabels
    axColors.set_yticklabels(legendLabels, fontdict=legendFont)
    
    #Fit layout for smaller image
    plt.tight_layout()
    
    #Save plot
    plt.savefig(outName,dpi=outRes)

if __name__=='__main__':
    
    #parameters for the function
    mapSize = (20,20)
    mapTitle = 'Cerrado Crop Health Index Comparison to 11y History (2006-2016) Arabica, '
    mapFile = ('/media/olivier/olivier_ext1/gedata_current/jde_coffee/MODIS/collection6/terra/Brazil/CER/'+
                'ndvi_2018-03-06_CompareToDecile_0BelowMin_110AboveMax_arabica_maskedbelow15pct.tif')
    backgroundLabel = 'Less than 15% Arabica'
    boundaryFile = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/data/Brazil/CER/aoi/AOI_CER_microregions'
    citiesFile = '/media/olivier/olivier_ext1/gedata_current/jde_coffee/data/Brazil/CER/places/CER_cities'
    citiesField = 'name'
    citiesLabelSize = 18
    citiesMarkerSize = 6
    notePosition = (0.05, 0.02)
    legendPosition = (0.05, 0.2)
    outName = '/home/olivier/Desktop/test.png'
    outRes = 200
    
    
    mapModisRanking(mapSize=mapSize, mapTitle=mapTitle, mapFile=mapFile, 
                    boundaryFile=boundaryFile, notePosition=notePosition, 
                    legendPosition=legendPosition, outName=outName, outRes=outRes,
                    backgroundLabel=backgroundLabel,
                    citiesFile=citiesFile, citiesField=citiesField,
                    citiesLabelSize=citiesLabelSize, citiesMarkerSize=citiesMarkerSize)
    