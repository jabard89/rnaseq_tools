"""
Created on 03/10/2023
@author: Jared Bard
fixes Skud gtf file downloaded from genbank to make it compatible with GenomicFeatures package in R
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import gffutils

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF3 file")
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export gff3 to")
	#args = parser.parse_args()
	args = parser.parse_args(["/home/jabard89/midway3-scratch/SK1_SK12/index/"
			                    "GCA_947243785.1_Skud-ZP591_genomic.gtf",                          
	                          	"/home/jabard89/midway3-scratch/SK1_SK12/index/"
			                    "GCA_947243785.1_Skud-ZP591_genomic.gtf"])
	
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
		
	db = gffutils.create_db(args.gff_file, ":memory:")

	with open(args.output_file, 'a') as fout:
		for d in db.directives:
			fout.write('## {0}\n'.format(d))

		for feature in db.all_features():
			feature.chrom = 'pombe'+feature.chrom
			fout.write(str(feature) + '\n')
		

	