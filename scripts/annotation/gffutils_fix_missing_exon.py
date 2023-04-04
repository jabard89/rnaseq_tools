"""
Created on 3/12/2023
@author: Jared Bard
Some genbank gtfs don't have exon features, which are needed by the R function makeTxDbFromGFF
Adding exons by copying CDS annotations, then adding a stop codon to the exon (required by makeTxDbFromGFF)
Some genes are also annotated as being in two parts:
gffutils can recognize these using the merge_strategy="merge" option
This function requires stop codons in the gtf file as currently written to identify the last exon
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
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/Data_JB/RNA-seq/SK1_SK12/GCA_947243785.1_Skud-ZP591_genomic.gtf",                          
	#                           	"/home/jabard89/Dropbox/Data_JB/RNA-seq/SK1_SK12/GCA_947243785.1_Skud-ZP591_genomic.cleaned.gtf"])
	
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
	
	# need to use merge here because the tRNA features are divided into multiple listings by "part"
	# but with the same gene_ids, which the import function does not like at all
	# also need to use gtf_subfeature='cds' because gtf doesn't have transcript features or exon features
	# this should reconstruct them
	db = gffutils.create_db(args.gff_file, ":memory:",id_spec={'gene':'gene_id','transcript':'Transcript'},
		merge_strategy="merge",gtf_subfeature='CDS')

	with open(args.output_file, 'a') as fout:
		for d in db.directives:
			fout.write('## {0}\n'.format(d))

		for feature in db.all_features():
			if feature.featuretype == 'CDS':
				fout.write(str(feature) + '\n') #first write the original CDS feature before changing into an exon
				# exons include the stop codons, but CDS do not. Need to add 3 nt for stop codon to the final exon
				# first grab parent gene and use it to find stop codon
				parent_gene = db[feature['gene_id'][0]]
				try:
					stop = list(db.children(parent_gene,featuretype='stop_codon'))[0]
				except IndexError:
					print("# Error: no stop codon found for gene {}".format(parent_gene['gene_id'][0]))
					sys.exit(1)
				if feature.strand == "+":
					if feature.end == stop.end-3:
						feature.end =feature.end+3
				else:
					if feature.start == stop.start+3:
						feature.start=feature.start-3
				feature.featuretype='exon'
				fout.write(str(feature) + '\n')
			else:
				fout.write(str(feature) + '\n')
		

	