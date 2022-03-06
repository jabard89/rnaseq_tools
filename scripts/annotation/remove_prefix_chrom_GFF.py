"""
Created on 9/27/2021
Edited 3/4/2022 to remove prefixes
Removes chrom prefix (chrI -> I). Keeps the original comments at the tops and adds new ones
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
	parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF sequence")
	parser.add_argument(dest="chrom_prefix",default="chr",type=str,help="prefix to add before chromosome")
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export GTF to")
	args = parser.parse_args()

	#Read in comments from original file
	header = []
	with open(args.gff_file,'r') as f:
		for line in f:
			if line[0].isalpha() or line[0].isdigit():
				break
			header.append(line)

	# Write out parameters
	with open(args.output_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))

		for line in header:
			f.write(line)

	# rename the chromosomes!
	gtf_file = HTSeq.GFF_Reader(args.gff_file,end_included=True)
	gtf_out_rename = []
	for i,feature in enumerate(gtf_file):
		chrom = feature.iv.chrom
		feature.iv.chrom = chrom.split(args.chrom_prefix,1)[1]
		gtf_out_rename.append(feature.get_gff_line())
	with open(args.output_file,'a') as f:
		f.writelines(gtf_out_rename)
