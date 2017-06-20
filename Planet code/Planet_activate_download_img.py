import requests
import shutil

# URL
url = "https://api.planet.com/data/v1/item-types/"

#API key
key = '77642b14ec2e4618a5688b47ea8db458'
#API user
user = 'olivierpfrancois@ge-data.com'
#API pwd
pwd = 'Gedata2015'

#Image information
img_link = ['https://api.planet.com/data/v1/item-types/REOrthoTile/items/20160821_133742_2329927_RapidEye-3/assets/',
	'https://api.planet.com/data/v1/item-types/PSOrthoTile/items/264120_2429901_2016-10-09_0e3a/assets/',
	'https://api.planet.com/data/v1/item-types/PSOrthoTile/items/221462_2429901_2016-08-08_0c41/assets/',
	'https://api.planet.com/data/v1/item-types/PSOrthoTile/items/264120_2329927_2016-10-09_0e3a/assets/',
	'https://api.planet.com/data/v1/item-types/PSOrthoTile/items/221462_2329827_2016-08-08_0c41/assets/',
	'https://api.planet.com/data/v1/item-types/PSOrthoTile/items/221462_2329928_2016-08-08_0c41/assets/',
	'https://api.planet.com/data/v1/item-types/PSOrthoTile/items/264124_2329728_2016-10-09_0e0d/assets/']
img_link = ['https://api.planet.com/data/v1/item-types/PSOrthoTile/items/264124_2329728_2016-10-09_0e0d/assets/']
img_id = ['20160821_133742_2329927_RapidEye-3',
	'264120_2429901_2016-10-09_0e3a',
	'221462_2429901_2016-08-08_0c41',
	'264120_2329927_2016-10-09_0e3a',
	'221462_2329827_2016-08-08_0c41',
	'221462_2329928_2016-08-08_0c41',
	'264124_2329728_2016-10-09_0e0d']
img_id = ['264124_2329728_2016-10-09_0e0d']
# Type of image: "REOrthoTile" (RapidEye), "PSOrthoTile" (PlanetScope)
#img_type = "PSOrthoTile"
#Format of image to donwload (available as raw image or simply for visuatlization)
#Can be visual or analytic
asset_type = "analytic"

# setup auth
session = requests.Session()
#session.auth = (user, pwd)
session.auth = (key, '')

for l in range(len(img_link)):
	# request an item
	item = \
		session.get(img_link[l])
		#session.get((url +"{}/items/{}/assets/").format(img_type, img_id))
	
	if item.json()[asset_type]['status']=='inactive':
		print 'item '+str(l)+': '+img_id[l]+' is inactive, sending activation request'
		# extract the activation url from the item for the desired asset
		img_activation_url = item.json()[asset_type]["_links"]["activate"]
		
		# request activation
		response = session.post(img_activation_url)
		print 'activation response status: '+str(response.status_code)
		
		item = \
			session.get(img_link[l])
		
		if item.json()[asset_type]['status']=='active':
			response = requests.get(item.json()[asset_type]['location'], stream=True)
			with open(img_id[l]+'.tif', 'wb') as out_file:
				shutil.copyfileobj(response.raw, out_file)
	
	if item.json()[asset_type]['status']=='activating':
		print 'item '+str(l)+': '+img_id[l]+' is activating'
	
	if item.json()[asset_type]['status']=='active':
		print 'item '+str(l)+': '+img_id[l]+' is active, starting download'
		response = requests.get(item.json()[asset_type]['location'], stream=True)
		with open(img_id[l]+'.tif', 'wb') as out_file:
			shutil.copyfileobj(response.raw, out_file)

#print item.json()[asset_type]
	
del response
