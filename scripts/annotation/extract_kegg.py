# Interactive script for accessing kegg pathways
# Jared Bard 07/13/2022
# Get organism code from https://www.kegg.jp/kegg/catalog/org_list.html
# API info https://www.kegg.jp/kegg/rest/keggapi.html
import pandas as pd

def retrieve_brite_dict(brite_id):
	"""
	Given a brite number, this will return a flat dictionary with every gene in every category
	"""
	import requests
	BASE = 'https://rest.kegg.jp'
	OPERATION = 'get'
	URL = BASE + '/' + OPERATION + '/br:' + brite_id + '/json'
	results = requests.get(URL)
	if not results.ok:
		print('Something went wrong in downloading the entry',results.status_code)
	entry = results.json()
	print('Retreived ' + entry['name'])
	return(entry)

t = retrieve_brite_dict("hsa03012")
	
def get_all_levels(entry):
	# returns a level_dict where each entry contains the nested dictionary below that level
	categories = {}
	def find_levels(l,parent=None):
		if 'children' in l.keys():
			if parent:
				name = parent + '.' + l['name']
			else:
				name = l['name']
			categories[name] = l
			for c in l['children']:
				find_levels(c,name)
	find_levels(entry)
	return categories

t_flat = get_all_levels(t)

def get_all_genes(entry):
	# given a dictionary (any item in the dictionary returned by get_all_levels)
	# returns a dataframe containing all the genes in that dictionary
	# including the genes in children
	# check levels using level_dict.keys()

	out_tuples = [] # (ID,KO)
	def extract_entries(l):
		for c in l['children']:
			if 'children' in c.keys():
					extract_entries(c)
			else:
				temp = c['name']
				split = temp.split('\t')
				gene = split[0].split()[0]
				ko = split[1].split()[0]
				out_tuples.append((gene,ko))
	extract_entries(entry)
	return pd.DataFrame(out_tuples,columns=['ORF','KO'])

hsa_eIF1 = get_all_genes(t_flat['hsa03012.Eukaryotic type.Initiation factors.eIF-1'])
b = get_all_genes(a['sce03012.Eukaryotic type.Initiation factors.eIF-1'])
