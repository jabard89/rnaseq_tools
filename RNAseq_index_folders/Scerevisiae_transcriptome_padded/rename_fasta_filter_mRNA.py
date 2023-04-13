# Written by: Jared Bard
# Date: 4/13/2023
# uses Biopython to rename fasta sequences to just the ORF name and filter for mRNA

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="in_file",default=None,type=str,help="input fasta file")
	parser.add_argument(dest="out_file",default=None,type=str,help="output fasta file")
	#args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/RNAseq_index_folders/"
	# 						"Scerevisiae_transcriptome_padded/Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_30down.fasta",""])
	args = parser.parse_args()

	#Read in comments from original file
	header = []
	with open(args.in_file,'r') as f:
		for line in f:
			if line[0] != "#":
				break
			header.append(line)

	# Write out parameters
	with open(args.out_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))

		for line in header:
			f.write(line)

	# input fastas
	input_fasta=list(SeqIO.parse(args.in_file,format='fasta'))
	output_fasta=[record for record in input_fasta]
	# rename sequences and filter for mRNA
	for record in output_fasta:
		try:
			ORF = re.search('transcript:(.+?)_mRNA', record.id).group(1)
		except AttributeError:
			print("Skipping {}".format(record.id))
			record=None
			continue
		record.id = ORF
		record.name =ORF
	with open(args.out_file,'a') as f:
		SeqIO.write(output_fasta,f,"fasta")