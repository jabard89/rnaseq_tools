# Written by: Jared Bard
# Date: 4/13/2023
# uses Biopython to rename fasta sequences to just the ORF name and filter for mRNA

import sys, os, math, random, argparse, csv, re
from Bio import SeqIO
from Bio import SeqRecord

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="in_file",default=None,type=str,help="input fasta file")
	parser.add_argument(dest="out_file",default=None,type=str,help="output fasta file")
	args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/RNAseq_index_folders/"
							"Scerevisiae_transcriptome_padded/Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_30down.fasta",""])
	args = parser.parse_args()

	# input fastas
	input_fasta=list(SeqIO.parse(args.in_file,format='fasta'))
	output_fasta=[]
	# rename sequences and filter for mRNA
	for record in input_fasta:
		try:
			ORF = re.search('transcript:(.+?)_mRNA', record.id).group(1)
		except AttributeError:
			print("Skipping {}".format(record.id))
			continue
		out_record = SeqRecord.SeqRecord(seq=record.seq,id=ORF,name=ORF,description=record.description)
		output_fasta.append(out_record)
	with open(args.out_file,'w') as f:
		SeqIO.write(output_fasta,f,"fasta")
