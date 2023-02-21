#! python
# Counts how many mismatches are in each gene, as output by hisat-3n_table
# Jared Bard
# 02/20/2023
# based on pile_reads_from_bam_genome_paired.py
# uses htseq to assign positions to genes, then goes through every line in the hisat-3n table and counts

import sys, os, math, random, argparse, csv, warnings, gzip, pathlib, re, itertools
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument('-w','--window',default=0,help="how many bp to use as window outside CDS")
	parser.add_argument(dest="in_gtf",type=pathlib.Path,default=None,help="input gtf file")
	parser.add_argument(dest="in_table", type=pathlib.Path,default=None, help="input conversion table file")
	parser.add_argument(dest="out_counts",type=pathlib.Path,default=None,help="name of the gene count output tsv")
	# Optional arguments
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/HeatShockQmRNA/JB_analysis/annotation/Saccharomyces_cerevisiae.R64-1-1.105_Spombe_geneid.gff3",
	# 						"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/230128_4TU/snake/tables/JB162_head.tsv",
	# 						"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/230128_4TU/snake/JB162_head_mismatch_pergene.tsv"])

	# Read input
	if not os.path.isfile(args.in_gtf):
		raise IOError("# Error: file {} does not exist".format(args.in_gtf))
	if not os.path.isfile(args.in_table):
		raise IOError("# Error: file {} does not exist".format(args.in_table))
	
	# check if either output is a compressed gz

	if os.path.isfile(args.out_counts):
		raise IOError("# Error: file {} already exists, please remove it".format(args.out_counts))
		
	with open(args.out_counts,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))

	gtf = HTSeq.GFF_Reader(str(args.in_gtf))
	
	# make a dictionary of genes with intervals that include the 5' and 3' UTRs
	print("Making database using {}".format(args.in_gtf))
	print("Exporting to {}".format(args.out_counts))
	gene_dict = {}
	for feature in gtf: # use a gtf where every item has a gene_id tag
		if feature.type in ['gene','CDS','5UTR','3UTR']:
			if 'gene_id' not in feature.attr:
				raise IOError("# Error: gtf {} should have a gene_id tag for every entry".format(args.in_gtf))
			name = feature.attr['gene_id']
		else:
			continue
		if name not in gene_dict:
			gene_dict[name] = {'gene':[],'5UTR':[],'CDS':[],'3UTR':[],'intron':[]}
		gene_dict[name][feature.type].append(feature)
	
	# remove genes with duplicate entries
	remove_list = []
	for gene in gene_dict:
		if len(gene_dict[gene]['gene']) != 1:
			remove_list.append(gene)
	for gene in remove_list:
		print('Removing {} because of duplicate entries or mismatched names'.format(gene))
		gene_dict.pop(gene)
		gene_ivs.pop(gene)
	
	# make a feature array with the intervals of all the CDS
	cds_array = HTSeq.GenomicArrayOfSets("auto",stranded=True)
	counts = {} # dictionary to store counts in
	for gene in gene_dict:
		counts[gene] = {'convertedBaseCount':0,'unconvertedBaseCount':0}
		if gene_dict[gene]['CDS']:
			for region in gene_dict[gene]['CDS']:
				cds_array[region.iv] += gene
		
	print("Reading {}".format(args.in_table))
	conversion_table = pd.read_csv(args.in_table,sep='\t')
	print("Counting table:")
	for index, row in conversion_table.iterrows():
		if index % 100000 == 0:
			print("{}/{}".format(index,len(conversion_table.index)))
		pos = HTSeq.GenomicPosition(chrom=row['ref'],pos=row['pos']-1,strand=row['strand'])
		set = cds_array[pos] #check if any genes are at the location
		if len(set)==1:
			counts[list(set)[0]]['convertedBaseCount']+=row['convertedBaseCount']
			counts[list(set)[0]]['unconvertedBaseCount']+=row['unconvertedBaseCount']
	
	counts_out = []
	for gene in counts:
		df = pd.DataFrame(counts[gene],index=[gene])
		df['ORF'] = df.index
		counts_out.append(df)
	pd.concat(counts_out).to_csv(args.out_counts,sep='\t',index=False,mode='a')
