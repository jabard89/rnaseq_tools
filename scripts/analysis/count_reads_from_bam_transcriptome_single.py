#! python
# Extracts reads mapped to every gene from data that was mapped to genomes
# Jared Bard
# 06/15/2022
# Modified version of pile_reads_from_bam_genome_single, but for counting reads mapped with STAR to a transcriptome fasta
# use for analyzing ribosome profiling data
# align reads 17 bp from the 5' end of the read
# can output either tsv or tsv.gz files

import sys, os, math, random, argparse, csv, warnings, gzip, pathlib, re, itertools
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq
from Bio import SeqIO


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument('-i',"--invert_strands",action="store_true",help="invert the strands when looking for matches")
	parser.add_argument(dest="in_fasta",type=pathlib.Path,default=None,help="input fasta file used for mapping")
	parser.add_argument(dest="in_bam", type=pathlib.Path,default=None,help="input bam file")
	parser.add_argument(dest="out_dist",type=pathlib.Path,default=None,help="name of the gene distribution output tsv")
	parser.add_argument(dest="out_counts",type=pathlib.Path,default=None,help="name of the gene count output tsv")
	parser.add_argument("-f","--trim_front",dest="trim_front",default=200,type=int,help="nt to trim from front")
	parser.add_argument("-b","--trim_back",dest="trim_back",default=30,type=int,help="nt to trim from back")
	parser.add_argument("-m","--min_length",dest="min_length",default=100,type=int,help="minimum remaining nt")
	# Optional arguments
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/hs-datasets/analysis/Jared/compare_riboseq/sgd_20210422_orf_coding_all.fasta",
	# 						"/home/jabard89/Dropbox/code_JB/repos/hs-datasets/analysis/Jared/compare_riboseq/Weinberg.2016_mrna_30C.control_BY4741_ribozero_BR1_transcriptome_Aligned.out.bam",
	# 						"/home/jabard89/Dropbox/code_JB/repos/hs-datasets/analysis/Jared/compare_riboseq/test_dist.tsv.gz",
	# 						"/home/jabard89/Dropbox/code_JB/repos/hs-datasets/analysis/Jared/compare_riboseq/test_counts.tsv"])

	# Read input
	if not os.path.isfile(args.in_fasta):
		raise IOError("# Error: file {} does not exist".format(args.in_fasta))
	if not os.path.isfile(args.in_bam):
		raise IOError("# Error: file {} does not exist".format(args.in_bam))
	
	# check if either output is a compressed gz
	dist_gzip_flag = False
	counts_gzip_flag = False
	dist_uncompressed = args.out_dist
	dist_compressed = None
	counts_uncompressed = args.out_counts
	counts_compressed = None
	if args.out_dist.suffix == ".gz":
		dist_gzip_flag = True
		dist_uncompressed = args.out_dist.parent / args.out_dist.stem
		dist_compressed = args.out_dist
	if args.out_counts.suffix == ".gz":
		counts_gzip_flag = True
		counts_uncompressed = args.out_counts.parent / args.out_counts.stem
		counts_compressed = args.out_counts

	for file in [dist_uncompressed,dist_compressed,counts_uncompressed,counts_compressed]:
		if file and os.path.isfile(file):
			raise IOError("# Error: file {} already exists, please remove it".format(file))
		
	with open(dist_uncompressed,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))
	with open(counts_uncompressed,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))	

	bam = HTSeq.BAM_Reader(str(args.in_bam))
	input_fasta=list(SeqIO.parse(args.in_fasta,format='fasta'))

	# make a dictionary of genes with intervals that include the 5' and 3' UTRs
	print("Making Databases for {} using {}".format(args.in_bam,args.in_fasta))
	print("Exporting to {} and {}".format(args.out_dist,args.out_counts))
	gene_names = [seq.name for seq in input_fasta]
	gene_lengths = [len(seq) for seq in input_fasta]
	gene_dict = {input_fasta[i].name: len(input_fasta[i]) for i in range(len(input_fasta))}

	dist = {}
	for i,gene in enumerate(gene_names):
		eff_len = gene_lengths[i] - args.trim_front - args.trim_back
		if eff_len >= args.min_length:
			transcript_pos = list(range(0,gene_lengths[i]))
			dist[gene] = {'index':transcript_pos,'Count':dict.fromkeys(transcript_pos,0),'eff_len':eff_len}


	counter = 0
	print("Counting reads:")
	for read in bam:
		counter += 1
		if counter % 100000 == 0:
			print(counter)
		if not read.aligned or read.optional_field("NH") > 1: # skip multi-mapping reads
			continue
		if args.invert_strands:
			read.iv = invert_strand(read.iv)
		if not read.iv.chrom in dist:
			continue
		if read.iv.start < args.trim_front or read.iv.start >= (gene_dict[read.iv.chrom]-args.trim_back):
			continue
		# assign the read to the correct position on the gene
		
		dist[read.iv.chrom]['Count'][read.iv.start] += 1
		
	print("Exporting counts!")
	# reformat distribution counts
	dist_out = []
	for gene in dist:
		df = pd.DataFrame.from_dict(dist[gene]['Count'],orient='index',columns=['Count'])
		df['Pos'] = dist[gene]['index']
		df = df[df.Count != 0]
		df['ORF'] = gene
		dist_out.append(df)
	pd.concat(dist_out).to_csv(dist_uncompressed,sep='\t',mode='a',index=False)
	if dist_gzip_flag:
		with open(dist_uncompressed,'rb') as src, gzip.open(dist_compressed,'wb') as dst:
			dst.writelines(src)
		os.remove(dist_uncompressed)
	
	counts_out = []
	for gene in dist:
		df = pd.DataFrame({'Count':sum(dist[gene]['Count'].values())},index=[gene])
		#df['length_trans']=len(dist[gene]['index'])
		df['CDS_length']=gene_dict[gene]
		df['Eff_length']=dist[gene]['eff_len']
		df['ORF'] = gene
		counts_out.append(df)
	counts_df = pd.concat(counts_out)
	counts_df['RPKB'] = (counts_df['Count'] / (counts_df['Eff_length']/1000))
	totalRPKB = counts_df['RPKB'].sum()
	counts_df['TPM'] = counts_df['RPKB']*1e6/totalRPKB
	counts_df.drop('RPKB',axis=1).to_csv(counts_uncompressed,sep='\t',index=False,mode='a')
	if counts_gzip_flag:
		with open(counts_uncompressed,'rb') as src, gzip.open(counts_compressed,'wb') as dst:
			dst.writelines(src)
		os.remove(counts_uncompressed)
