# Customize this starter script by adding code
# to the run_script function. See the Help for
# complete information on how to create a script
# and use Script Runner.

""" Your Description of the script goes here """

# Imports

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
import qgis.utils as qu
from datetime import datetime
import os


def run_script(iface):
    
    #########################################################################################################################
    #########################################################################################################################
    # #PARAMETERS
    
    # Satellite
    satelliteModis = 'terra'  # 'terra' # 'aqua'
    
    # MODIS date to map out
    modisDates = ['2017-12-03','2018-01-01','2018-02-02']
    
    # % of coffee masks to map out (the rasters for each should have been prepared in advance
    modisPct = ['5', '15']
    
    # Address for the working directory for MODIS and for the data (boundaries, cities)
    root = '/media/olivier/olivier_ext1/gedata_current/jde_coffee' # '/home/olivierp/jde_coffee'
    modisPrefix = root+'/MODIS/collection6/'+satelliteModis+'/Vietnam' 
    dataPrefix = root+'/data/Vietnam' 
        
    # Destination folder for the maps from the modisPrefix
    destFolder = 'maps'
    
    # Names of the States/regions to map out. There should be one entry per coffee variety if the state contains more than one
    states = ['LD']
    # Coffee varieties for each of the states
    varieties = ['coffee']  # Coffee variety for each map
    # Titles for each of the modis maps, to which the modis date will be added at the end
    mapTitles = ['Lam Dong Crop Health Index \nComparison to 11y History (2006-2016) \nCoffee, ']
    # Name of the boundary shapefile for each of the states
    boundaries = ['aoi/AOI_LD.shp']
    # Name of the shapefile with the cities for each of the states
    cities = ['LD_cities.shp']
    # Size of the maps for each of the states (options are 'PORTRAIT', 'LANDSCAPE' AND 'SQUARE')
    mapSizes = ['LANDSCAPE']
    # ##MAP BOXES Parameters
    # Parameters for the map boxes of each region
    topX = [12]  # x position (cm) of the top left corner of the map box
    topY = [12]  # y position (cm) of the top left corner of the map box
    # Size of the bottom margin below the map box. It will condition the size of the mapbox together with the top y position
    marginBottom = [10]
    # Should the map box be framed?
    frameMapBox = False
    
    # ##TITLE Parameters
    titleX = [-35]  # x position (cm) of the top left corner of the title box
    titleY = [4]  # y position (cm) of the top left corner of the title box
    titleHeight = [25]  # Height (cm) of the title box
    
    # ##LEGEND Parameters
    legendX = [6]  # x position (cm) of the top left corner of the legend box
    legendY = [25]  # y position (cm) of the top left corner of the legend box
    
    # ##COMMENT BOX Parameters 
    # The comment box is set to have the same width as the legend and to be set a little below the legend
    commentGap = [63]  # How much below the legend should the comment box start?
    commentHeight = 40  # Height of the comment box
    # Comment to add to each of the maps 
    comment = 'Notes:'\
                '\n--Most pixels contain other land uses besides coffee.'\
                '\n--The index shows health of vegetation in each pixel compared to reference years, '\
                'not coffee production directly. '\
                '\n--Vegetative health is affected mostly by '\
                'natural factors (rain, etc.) but can also be affected by human intervention (pruning, etc.).'
    
    # ##SCALE Parameters
    scaleX = 20  # How much to the left of the left side of the paper should the scale bar end
    scaleY = 15  # How much above the bottom of the paper should the scale be
    
    # Format for export
    exportFormat = 'PNG'  # 'PDF'
    # Resolution for export
    exportResolution = 300
    
    # Parameters transparency for the title and legend frames. Currently not used as frames are set with no background
    itemsTransparency1 = 190
    itemsTransparency2 = 100
    
    #########################################################################################################################
    #########################################################################################################################
    # #ACTIVE CODE
    
    # Assign map canvas to a variable
    canvas = qu.iface.mapCanvas()
    # Change map units to meters
    canvas.setMapUnits(QGis.Meters)
    
    # Get the registry
    registry = QgsMapLayerRegistry.instance()
    
    for s in range(len(states)):
        for modisDate in modisDates:
            for p in modisPct:
                
                #####################################################
                # LOAD LAYERS AND SET THE SYMBOLOGY
                
                # Declare the name for the boundary layer to use by QGIS
                boundNm = 'Less than ' + p + '% ' + varieties[s]
                
                # Remove previous layers (boundary, cities, modis) if they exist already
                for layer in registry.mapLayers().values():
                    if layer.name() in ['Cities', ''] or 'Less than ' in layer.name():
                        registry.removeMapLayers([layer.id()])
                
                # Import the shapefile with the cities
                city = QgsVectorLayer(dataPrefix + '/' + states[s] + '/places/' + cities[s], 'Cities', 'ogr')
                # Load style for the cities
                city.loadNamedStyle(root+'/MODIS/decile_comparison_style_cities.qml')
                
                # Import the modis raster layer
                modis = QgsRasterLayer(modisPrefix + '/' + states[s] + '/ndvi_' + modisDate + 
                                       '_CompareToDecile_0BelowMin_110AboveMax_' + varieties[s] + '_maskedbelow' + p + 'pct.tif')
                
                # Load style for modis
                modis.loadNamedStyle(root+'/MODIS/decile_comparison_style_purplebluescale.qml')
                
                # Import the shapefile with the state boundaries
                bound = QgsVectorLayer(dataPrefix + '/' + states[s] + '/' + boundaries[s], boundNm, 'ogr')
                # Change the style for the boundary layer
                # Prepare symbol properties
                properties = {}
                properties["color"] = '#000000'
                properties["color_border"] = '#000000'
                properties["width_border"] = '0.26'
                properties["style"] = 'no'
                # Create a new symbol with these properties
                sym = QgsFillSymbolV2.createSimple(properties)
                # Apply the symbol to the boundary layer
                bound.rendererV2().setSymbol(sym)
                # Repaint the layer to force an update
                # bound.triggerRepaint()
                
                # Add the layers to the registry
                registry.addMapLayers([modis, bound, city])
                
                # Zoom to extent of boundary layer
                canvas.setExtent(bound.extent())
                
                #####################################################
                # CREATE THE MAP
                
                # Set the size of the map
                if mapSizes[s] == 'PORTRAIT': 
                    paperSize = QSize(210, 297)     
                elif mapSizes[s] == 'LANDSCAPE': 
                    paperSize = QSize(297, 210)
                elif mapSizes[s] == 'SQUARE': 
                    paperSize = QSize(244, 244)
                # Set the size of the map box where the MODIS will show
                rendererSize = QSize(paperSize.width() - 2 * topX[s], paperSize.height() - topY[s] - marginBottom[s])
                composerMap_w = rendererSize.width()
                composerMap_h = rendererSize.height()
                
                # Initialize Map renderer
                mapRenderer = canvas.mapRenderer()
                mapRenderer.setOutputUnits(QgsMapRenderer.Millimeters)
                mapRenderer.setOutputSize(rendererSize, exportResolution)
                
                # Add layers to map renderer (the last added layer is the one on the bottom)
                # Get the layer id and names into a dictionary
                layers = {}
                for layer in registry.mapLayers().values():
                    layers[layer.name()] = layer.id()
                # Add them in the right order to the map
                mapRenderer.setLayerSet([layers[boundNm], layers['Cities'], layers['']])
                
                # Set the zoom for the map
                curExtent = modis.extent()
                mapRenderer.setExtent(curExtent)
                
                # Make Composition
                c = QgsComposition(mapRenderer)
                c.setPlotStyle(QgsComposition.Print)
                c.setPaperSize(paperSize.width(), paperSize.height())
                c.setPrintResolution(exportResolution)
                
                # Position the Map and add to Composition
                composerMap = QgsComposerMap(c, topX[s] , topY[s], composerMap_w, composerMap_h)
                composerMap.setFrameEnabled(frameMapBox)
                composerMap.setFrameOutlineWidth(0.5)  # Used only if a frame is enabled for the map box
                c.addItem(composerMap)
                
                #####################################################
                # ADD THE TITLE
                
                composerLabel = QgsComposerLabel(c)
                composerLabel.setFont(QFont('Ubuntu', 18))
                composerLabel.setFontColor(QColor(0, 0, 0))
                composerLabel.setHAlign(0x0004)  # # check Qt.AlignmentFlag in http://pyqt.sourceforge.net/Docs/PyQt4/qt.html#AlignmentFlag-enum
                # Size of title frame and placement: (startX, startY [from top], width, height)
                composerLabel.setSceneRect(QRectF(titleX[s], titleY[s], paperSize.width(), titleHeight[s]))
                date = datetime.strptime(modisDate, '%Y-%m-%d').date()  # Transform into date to be able to reformat
                composerLabel.setText(mapTitles[s] + date.strftime('%B %d, %Y'))
                # composerLabel.adjustSizeToText() #Should the box be adjusted to the text?
                composerLabel.setFrameEnabled(False)  # Should the title box be framed?
                composerLabel.setBackgroundEnabled(False)  # Should the title box have a background?
                composerLabel.setBackgroundColor(QColor(255, 255, 255, itemsTransparency1))  # Background color, used only if background is enabled
                c.addItem(composerLabel)
                
                #####################################################
                # ADD THE LEGEND
                
                # Initialize list of layers to show in legend 
                layerGroup = QgsLayerTreeGroup()
                layerGroup.addLayer(modis)
                layerGroup.addLayer(bound)
                
                legend = QgsComposerLegend(c)
                legend.setTitle("Legend")
                legend.modelV2().setRootGroup(layerGroup)
                legend.setStyleFont(QgsComposerLegendStyle.Title, QFont('Ubuntu', 16))
                legend.setStyleFont(QgsComposerLegendStyle.Subgroup, QFont('Ubuntu', 6))
                legend.setStyleFont(QgsComposerLegendStyle.SymbolLabel, QFont('Ubuntu', 12))
                legend.setBoxSpace(2)
                legend.setFrameEnabled(False)
                legend.setFrameOutlineWidth(0.25)  # Frame line size, used only if a frame is enabled
                legend.setFrameOutlineColor(QColor(0, 0, 0))  # Frame line color, used only if a frame is enabled
                legend.setBackgroundEnabled(False)
                legend.setBackgroundColor(QColor(255, 255, 255, itemsTransparency1))  # Background color, used only if background is enabled
                legend.setSymbolHeight(6)
                legend.setResizeToContents(True)        
                legendSize = legend.paintAndDetermineSize(None)
                # Set the position and size of the legend
                legend.setItemPosition(legendX[s], legendY[s], legendSize.width(), legendSize.height(), QgsComposerItem.UpperLeft, False, -1)
                c.addItem(legend)
                
                #####################################################
                # ADD THE COMMENT
                
                composerLabel = QgsComposerLabel(c)
                composerLabel.setFont(QFont('Ubuntu', 9))
                composerLabel.setFontColor(QColor(0, 0, 0))
                composerLabel.setHAlign(0x0001)  # # check Qt.AlignmentFlag in http://pyqt.sourceforge.net/Docs/PyQt4/qt.html#AlignmentFlag-enum
                # Size of comment frame and placement: (startX, startY [from top], width, height)
                composerLabel.setSceneRect(QRectF(legendX[s], legendY[s] + legendSize.height() + commentGap[s], legendSize.width(), commentHeight))
                composerLabel.setText(comment)
                
                composerLabel.setFrameEnabled(False)
                composerLabel.setBackgroundEnabled(False)
                composerLabel.setBackgroundColor(QColor(255, 255, 255, itemsTransparency1))
                c.addItem(composerLabel)
                
                #####################################################
                # ADD THE SCALE BAR
                
                scaleBar = QgsComposerScaleBar(c)
                scaleBar.setStyle('Double Box')
                scaleBar.setComposerMap(composerMap)
                scaleBar.setFrameEnabled(False)
                scaleBar.setBackgroundEnabled(False)
                scaleBar.setBackgroundColor(QColor(255, 255, 255, itemsTransparency2))
                scaleBar.setUnitLabeling('km')
                scaleBar.setNumMapUnitsPerScaleBarUnit(25)
                scaleBar.applyDefaultSize()
                scaleBar_x = paperSize.width() - ((scaleBar.segmentMillimeters() * scaleBar.numSegments())) - scaleX
                scaleBar_y = paperSize.height() - scaleY
                scaleBar.setItemPosition(scaleBar_x, scaleBar_y)
                scaleBar.setPen(QPen(QBrush(Qt.black), 0.2, Qt.SolidLine))
                scaleBar.setNumSegmentsLeft(0)
                scaleBar.setNumSegments(4)
                scaleBar.setFont(QFont('Ubuntu', 12))
                scaleBar.setHeight(3)
                scaleBar.setLabelBarSpace(1)
                c.addItem(scaleBar)
                
                '''
                # Final frame around map to hide artefacts of title and scalebar
                frameRect = QgsComposerShape(marginSize, marginSize, composerMap_w, composerMap_h, c)
                frameRect.setFrameOutlineWidth(0.5)
                frameRect.setFrameOutlineColor(QColor(0,0,0))
                frameRect.setBackgroundEnabled(False)
                frameRect.setShapeType(1)
                c.addItem(frameRect)
                '''
                
                #####################################################
                # EXPORT THE MAP
                
                # Create name
                outNm = (modisPrefix + '/' + destFolder + 
                         '/decile_comparison_' + states[s] + 
                         '_' + varieties[s] + '_' + p + 'pct_' + 
                         modisDate + '.' + exportFormat.lower())
                exportMap(c, exportFormat, outNm)
            
            # Remove raster layer (modis)
            # for layer in registry.mapLayers().values():
            #    if layer.name() == '':
            #        registry.removeMapLayers([layer.id()])
            
            # Clean the xml files
            xml = [f for f in os.listdir(os.path.join(modisPrefix, states[s])) if f.endswith('.aux.xml')]
            for f in xml:
                os.remove(os.path.join(modisPrefix + '/' + states[s], f))
        

def exportMap(theComposition, exportFormat, outFName):
    """
    ---------------------------------------------------------------------------------------------
    Function : exports the Map in PNG or PDF format and saves on disk
    ---------------------------------------------------------------------------------------------
    """
    # PDF Output
    if exportFormat == 'PDF':
        printer = QPrinter()
        printer.setOutputFormat(QPrinter.PdfFormat)
        printer.setOutputFileName(outFName)
        printer.setPaperSize(QSizeF(theComposition.paperWidth(), theComposition.paperHeight()), QPrinter.Millimeter)
        printer.setFullPage(True)
        printer.setColorMode(QPrinter.Color)
        printer.setResolution(theComposition.printResolution())

        pdfPainter = QPainter(printer)
        paperRectMM = printer.pageRect(QPrinter.Millimeter)
        paperRectPixel = printer.pageRect(QPrinter.DevicePixel)
        theComposition.render(pdfPainter, paperRectPixel, paperRectMM)
        pdfPainter.end()

    # IMAGE output
    elif exportFormat == 'PNG':

        # Set dimensions and resolution
        dpi = theComposition.printResolution()
        dpmm = (dpi / 25.4)
        width = int(dpmm * theComposition.paperWidth())
        height = int(dpmm * theComposition.paperHeight())

        # Create output image and initialize it
        image = QImage(QSize(width, height), QImage.Format_ARGB32)
        image.setDotsPerMeterX(dpmm * 1000)
        image.setDotsPerMeterY(dpmm * 1000)
        image.fill(0)

        # Render composition
        imagePainter = QPainter(image)
        sourceArea = QRectF(0, 0, theComposition.paperWidth(), theComposition.paperHeight())
        targetArea = QRectF(0, 0, width, height)
        theComposition.render(imagePainter, targetArea, sourceArea)
        imagePainter.end()

        # Save image to disk 
        image.save(outFName, exportFormat.lower())   
