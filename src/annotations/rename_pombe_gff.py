"""
Created on 11/23/2022
@author: Jared Bard
Add a gene_id to the pombe gff
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
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export GTF to")
	args = parser.parse_args()
	#args = parser.parse_args(["../annotation/Schizosaccharomyces_pombe_all_chromosomes.gff3",
	#                           "../annotation/Schizosaccharomyces_pombe_all_chromosomes_geneid.gff3"])

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
	gff_out = []
	# first scan through the gtf_file and pull all of the genes
	for feature in gff_file:
		if feature.type == "gene":
			feature.attr['gene_id'] = feature.name
			gff_out.append("###\n")
			gff_out.append(feature.get_gff_line())
		else:
			gff_out.append(feature.get_gff_line())

	# write the file!
	with open(args.output_file,'a') as f:
		f.writelines(gff_out)

