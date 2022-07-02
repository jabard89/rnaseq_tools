#! python
# Extracts gene lengths from gff including utrs if available
# optional, extract lengths of transcripts including the introns
# Jared Bard
# 06/16/2022
# adapated from pile_reads_from_bam_genome_single.py

import sys, os, math, random, argparse, csv, warnings, gzip, pathlib, re, itertools
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest="in_gtf",type=pathlib.Path,default=None,help="input gtf file")
	parser.add_argument(dest="out_file",type=pathlib.Path,default=None,help="name of the output tsv")
	# Optional arguments
	parser.add_argument("-i","--intron",dest="intron",action="store_true",help="include introns in gene lengths")
	args = parser.parse_args()
	
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/annotations/Scerevisiae.R64-1-1.104.yeastss.pelechano.gtf",
	# 						"/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/annotations/Scerevisiae.R64-1-1.104.yeastss.pelechano_lengths.tsv",
	# 						"-i"])
	# Read input
	if not os.path.isfile(args.in_gtf):
		raise IOError("# Error: file {} does not exist".format(args.in_gtf))
	
	with open(args.out_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))	

	gtf = HTSeq.GFF_Reader(str(args.in_gtf))
	
	# make a dictionary of genes with intervals that include the 5' and 3' UTRs
	gene_dict = {}
	gene_ivs = {}
	for feature in gtf:
		if feature.type == "gene":
			name = feature.attr['gene_id']
		elif feature.type in ['CDS','5UTR','3UTR']:
			name = feature.attr['protein_id']
		else:
			continue
		if name not in gene_dict:
			gene_dict[name] = {'gene':[],'5UTR':[],'CDS':[],'3UTR':[],'intron':[]}
			gene_ivs[name] = None
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

	#extract the total interval for each gene, including 3' and 5'
	for gene in gene_dict:
		tempStart = gene_dict[gene]['gene'][0].iv.start
		tempEnd = gene_dict[gene]['gene'][0].iv.end
		tempStrand = gene_dict[gene]['gene'][0].iv.strand
		tempChrom = gene_dict[gene]['gene'][0].iv.chrom
		for type in gene_dict[gene]:
			for item in gene_dict[gene][type]:
				if item.iv.start < tempStart:
					tempStart = item.iv.start
				if item.iv.end > tempEnd:
					tempEnd = item.iv.end	
		gene_ivs[gene] = HTSeq.GenomicInterval(tempChrom,tempStart,tempEnd,tempStrand)
	
	
	genes_out=[gene for gene in gene_ivs]
	if args.intron:
		lengths_out=[gene_ivs[gene].length for gene in genes_out]
	else:
		lengths_out=[]
		for gene in genes_out:
			strand = gene_dict[gene]['gene'][0].iv.strand
			transcript_parts = [] # parts of the transcript, 0=5UTR,1-n=CDS,n+1=3UTR
			if len(gene_dict[gene]['5UTR']) == 1:
				iv1_start = gene_dict[gene]['5UTR'][0].iv.start_d
				iv1_end = gene_dict[gene]['5UTR'][0].iv.end_d
			else:
				iv1_start = gene_ivs[gene].start_d
				iv1_end = gene_dict[gene]['gene'][0].iv.start_d
			if strand == '-':
				transcript_parts.append(range(iv1_start,iv1_end,-1))
			else:
				transcript_parts.append(range(iv1_start,iv1_end))

			for CDS in gene_dict[gene]['CDS']:
				if strand == '-':
					transcript_parts.append(range(CDS.iv.start_d,CDS.iv.end_d,-1))
				else:
					transcript_parts.append(range(CDS.iv.start_d,CDS.iv.end_d))
			
			if len(gene_dict[gene]['3UTR']) == 1:
				iv3_start = gene_dict[gene]['3UTR'][0].iv.start_d
				iv3_end = gene_dict[gene]['3UTR'][0].iv.end_d
			else:
				iv3_start = gene_dict[gene]['gene'][0].iv.end_d
				iv3_end = gene_ivs[gene].end_d
			if strand == '-':
				transcript_parts.append(range(iv3_start,iv3_end,-1))
			else:
				transcript_parts.append(range(iv3_start,iv3_end))
			transcript_len = sum([len(part) for part in transcript_parts]) #calculate the length of the CDS
			lengths_out.append(transcript_len)
	
	
	df_out = pd.DataFrame({'ORF':genes_out,'Length':lengths_out})
	df_out.to_csv(args.out_file,sep='\t',mode='a',index=False)
			