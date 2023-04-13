import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="in_file",default=None,type=str,help="input fasta file")
	parser.add_argument(dest="out_file",default=None,type=str,help="output fasta file")
	parser.add_argument("-f","--trim_front",dest="trim_front",default=200,type=int,help="nt to trim from front")
	parser.add_argument("-b","--trim_back",dest="trim_back",default=30,type=int,help="nt to trim from back")
	parser.add_argument("-m","--min_length",dest="min_length",default=100,type=int,help="minimum remaining nt")
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

	# input and trim fastas
	input_fasta=list(SeqIO.parse(args.in_file,format='fasta'))
	output_fasta=[record for record in input_fasta if len(record)>(args.trim_front+args.trim_back+args.min_length)]
	if args.trim_back == 0:
		output_fasta=[record[args.trim_front:] for record in output_fasta]
	else:
		output_fasta=[record[args.trim_front:-args.trim_back] for record in output_fasta]
	with open(args.out_file,'a') as f:
		SeqIO.write(output_fasta,f,"fasta")
