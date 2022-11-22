# script for pulling annotations from kegg and outputting them as list of ORFs
# Jared Bard 10/25/2022
# uses functions in extract_kegg

import pandas as pd
import os as os
import extract_kegg as ek


from extract_kegg import retrieve_brite_dict
from extract_kegg import get_all_levels
from extract_kegg import get_all_genes

ribi_brite_id = 'sce03009'
ribi_dict = ek.get_all_levels(ek.retrieve_brite_dict(ribi_brite_id))
a = ribi_dict['sce03009']
ribi_brite_list = ek.get_all_genes(a)