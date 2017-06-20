import os
import requests
from requests.auth import HTTPBasicAuth
from jq import jq
from shapely.geometry import mapping, Polygon, shape
import fiona
from fiona.crs import from_epsg

#Planet API key
key = '77642b14ec2e4618a5688b47ea8db458'
#Planet API user
user = 'olivierpfrancois@ge-data.com'
#Planet API pwd
pwd = 'Gedata2015'


# Max % cloud coverage on images
cloud = 0.3
# Type of image: "REOrthoTile" (RapidEye), "PSOrthoTile" (PlanetScope), 
# "PSScene3Band" (PlanetScope indivudal images)
img_type = "REOrthoTile"
# Date range
dates = {
'start': '2016-01-01',
'end': '2016-10-15'
}

# Method for coordinates of AOI
input_coord = 'shape' # 'shape' or 'manual'
# Provide the full address of the shapefile with the AOI
shp = '/media/olivier/olivier_ext/gedata_current/jde_coffee/data/MG/ZM/AOI/AOI_MG_ZM_geo.shp'
#'/home/olivier/Desktop/test_poly.shp'
# Provide coordinates if the method is manual
if input_coord=='manual':
	# Coordinates of box(es) to search: Should be at least 5 coordinates 
	# for each box. The first and last coordinates should be the same.
	# If more than one box, should be a list of lists
	coords = [
	[[-122.52227783203125,40.660847697284815],
	[-122.52227783203125,40.987154933797335],
	[-122.01690673828124,40.987154933797335],
	[-122.01690673828124,40.660847697284815],
	[-122.52227783203125,40.660847697284815]]
	]
elif input_coord=='shape':
	#Import shapefile
	coords = fiona.open(shp)
	coords = [p['geometry']['coordinates'] for p in coords]
	coords = [p[0] for p in coords] 

for box in range(len(coords)):
	
	# Create a json geometry object of box of interest
	geo_json_geometry = {
	"type": "Polygon",
	"coordinates": [coords[box]]
	}

	# filter for items the overlap with our chosen geometry
	geometry_filter = {
	"type": "GeometryFilter",
	"field_name": "geometry",
	"config": geo_json_geometry
	}
	
	# filter images acquired in a certain date range
	date_range_filter = {
	"type": "DateRangeFilter",
	"field_name": "acquired",
	"config": {
		"gte": dates['start']+"T00:00:00.000Z",
		"lte": dates['end']+"T00:00:00.000Z"
	}
	}
	
	# filter any images which are more than 50% clouds
	cloud_cover_filter = {
	"type": "RangeFilter",
	"field_name": "cloud_cover",
	"config": {
		"lte": cloud
	}
	}
	
	# create a filter that combines our geo and date filters
	# could also use an "OrFilter"
	images_filter = {
	"type": "AndFilter",
	"config": [geometry_filter, date_range_filter, cloud_cover_filter]
	}
	
	## use the filter to query the stats endpoint, this will give us a 
	# date bucketed histogram to show us how many items match our filter
	# Stats API request object
	#stats_endpoint_request = {
	#  "interval": "day",
	#  "item_types": [img],
	#  "filter": images_filter
	#}
	## fire off the POST request
	#result = \
	#  requests.post(
	#    'https://api.planet.com/data/v1/stats',
	#    auth=HTTPBasicAuth(user, pwd),
	#    json=stats_endpoint_request)
	
	
	# Search API request object
	search_endpoint_request = {
	"item_types": [img_type],
	"filter": images_filter
	}
	
	result = \
	requests.post(
		'https://api.planet.com/data/v1/quick-search',
		#auth=HTTPBasicAuth(user, pwd),
		auth=HTTPBasicAuth(key, ''),
		json=search_endpoint_request)
	
	#print result.text
	#print jq(".features[]").transform(result.json(), text_output=True)
	#Extract the IDs of the images
	ids = jq(".features[].id").transform(result.json(), text_output=True)
	ids = [str(i) for i in ids.replace('"','').split('\n')]
	_link = jq(".features[]._links.assets").transform(result.json(), text_output=True)
	_link = [str(i) for i in _link.replace('"','').split('\n')]
	cover = jq(".features[].properties.cloud_cover").transform(result.json(), text_output=True)
	cover = [str(i) for i in cover.replace('"','').split('\n')]
	reso = jq(".features[].properties.pixel_resolution").transform(result.json(), text_output=True)
	reso = [str(i) for i in reso.replace('"','').split('\n')]
	date = jq(".features[].properties.acquired").transform(result.json(), text_output=True)
	date = [str(i)[0:10] for i in date.replace('"','').split('\n')]
	azimuth = jq(".features[].properties.sun_azimuth").transform(result.json(), text_output=True)
	azimuth = [str(i) for i in azimuth.replace('"','').split('\n')]
	epsg = jq(".features[0].properties.epsg_code").transform(result.json())
	poly = jq(".features[].geometry").transform(result.json(), multiple_output=True)
	poly = [shape(l) for l in poly]
	
	#fmt = '{:<4}{:>35}{:>15}{:>15}{:>15}' #Declare the format for the output to be aligned
	#print(fmt.format('', 'Id', 'Cloud cover', 'Resolution', 'Sun azimuth'))
	#for j, (i, c, r, a) in enumerate(zip(ids, cover, reso, azimuth)):
	#    print(fmt.format(j, i, c, r, a))
	
	# Define a polygon feature geometry with the attributes
	planetschema = {
		'geometry': 'Polygon',
		'properties': {'id':'str', 'cover':'float', 'reso':'float', 'azimuth':'float', 'date':'str', 'link':'str'},
	}
	
	# Write a new Shapefile using fiona # crs=from_epsg(epsg),
	with fiona.open(img_type+'_'+dates['start']+'_'+dates['end']+'_'+'aoi'+str(box+1)+'.shp', 'w', crs=from_epsg(4326), driver='ESRI Shapefile', schema=planetschema) as c:
		for p in range(len(poly)):
			c.write({
				'geometry': mapping(poly[p]),
				'properties': {'id': ids[p], 'cover': cover[p], 'reso': reso[p], \
					'azimuth': azimuth[p], 'date': date[p], 'link': _link[p]},
			})
