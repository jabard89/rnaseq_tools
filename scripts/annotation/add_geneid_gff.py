"""
Created on 11/23/2022
@author: Jared Bard
For use with STAR --quantMode, add a gene_id attribute to every feature in a group
assumes that gff files can have two levels of nesting, but that features can have multiple parents!
strips quotes from gene_ids
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import HTSeq
import itertools
import re

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF3 file")
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export gff3 to")
	args = parser.parse_args()
	#args = parser.parse_args(["../src/Schizosaccharomyces_pombe_all_chromosomes.gff3",
	#                          "../src/Spombe_test.gff3"])

	for file in [args.gff_file]:
		if not os.path.isfile(file):
			raise IOError("# Error: file {} does not exist".format(file))	
	
	for file in [args.output_file]:
		if os.path.isfile(file):
			raise IOError("# Error: file {} already exists, please remove it".format(file))
	
	# Write out parameters
	with open(args.output_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))

	# read the GTF file, grouping together all the CDS by gene
	gff_file = HTSeq.GFF_Reader(args.gff_file,end_included=True)
	print_flag = False
	f0_dict = {}
	f1_dict = {} # dictionary of features that have genes as parents
	f2_dict = {} # dictionary of features with f1 parents
	renamed_features = {} # dictionary of features that were renamed

	# first scan through the gtf_file and pull all of the genes
	for feature in gff_file:
		if 'Parent' not in feature.attr:
			if 'gene_id' in feature.attr:
				gene_id = feature.attr['gene_id'].replace('"','')
			else:
				print("No gene_id found for {}, using ID".format(feature.name))
				gene_id = feature.name.replace('"','')
			# now add the gene to the gene dictionary
			f0_dict[feature.name] = {'gene_id':gene_id,'feature':feature,'children':[]}
	# scan for mRNAs associated with each gene
	for feature in gff_file:
		if ('Parent' in feature.attr) and (feature.attr['Parent'] in f0_dict):
			gene_id = f0_dict[feature.attr['Parent']]['gene_id']
			feature.attr['gene_id']=gene_id
			renamed_features[feature.name] = feature
			f0_dict[feature.attr['Parent']]['children'].append(feature.name)
			f1_dict[feature.name] = {'gene_id':gene_id,'feature':feature,'children':[]}
	# scan for other features
	for feature in gff_file:
		if ('Parent' in feature.attr) and (feature.attr['Parent'] in f1_dict):
			gene_id=f1_dict[feature.attr['Parent']]['gene_id']
			feature.attr['gene_id']=gene_id
			renamed_features[feature.name] = feature
			f1_dict[feature.attr['Parent']]['children'].append(feature.name)
			f2_dict[feature.name] = {'gene_id':gene_id,'feature':feature}
		else:
			if ('Parent' in feature.attr):
				if (feature.attr['Parent'] not in f0_dict) and (feature.attr['Parent'] not in f1_dict):
					print("Can't find parent for {}".format(feature.name))

	
	# now remake the gff
	gff_out = []
	for feature in gff_file:
		if feature.type == "gene":
			gff_out.append("###\n")
			gff_out.append(feature.get_gff_line())
		elif feature.name in renamed_features:
			gff_out.append(renamed_features[feature.name].get_gff_line())
		else:
			gff_out.append(feature.get_gff_line())

	# write the file!
	with open(args.output_file,'a') as f:
		f.writelines(gff_out)






