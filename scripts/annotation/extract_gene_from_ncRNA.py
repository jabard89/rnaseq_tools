# Jared Bard 220616
# take a fasta input from SGD and extract the first word of the description as the gene name
import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="in_file",default=None,type=str,help="input fasta file")
	parser.add_argument(dest="out_file",default=None,type=str,help="output tsv file")
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/genomes/sgd_20210427_rna_coding.fasta",
	# 						"/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/annotations/sgd_20210427_rna_genes.tsv"])
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

	# input and trim fastas
	input_fasta=list(SeqIO.parse(args.in_file,format='fasta'))
	ORFs_out = []
	names_out = []
	for item in input_fasta:
		ORFs_out.append(item.name)
		names_out.append(item.description.split()[1])
	df_out = pd.DataFrame({'ORF':ORFs_out,'gene':names_out})
	df_out.to_csv(args.out_file,sep='\t',mode='a',index=False)