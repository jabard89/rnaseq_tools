#! python
# Counts the number of reads mapping to ncRNAs from directional paired reads
# Jared Bard
# 06/19/2022
# based on pile_reads_from_bam_genome_paired.py

import sys, os, math, random, argparse, csv, warnings, gzip, pathlib, re, itertools
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument('-i','--invert_strands',action="store_true",help="invert the strands when looking for matches")
	parser.add_argument(dest="in_gtf",type=pathlib.Path,default=None,help="input gtf file")
	parser.add_argument(dest="in_bam", type=pathlib.Path,default=None, help="input bam file")
	parser.add_argument(dest="out_counts",type=pathlib.Path,default=None,help="name of the gene count output tsv")
	# Optional arguments
	# args = parser.parse_args()
	args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/annotations/Scerevisiae.R64-1-1.104.yeastss.pelechano.gtf",
							"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/220301/star/mapped_reads/220403_STAR_HG49_Aligned.sortedByCoord.dedup.out.bam",
							"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/220301/star/HG49/HG49_counts_test.tsv"])

	# Read input
	if not os.path.isfile(args.in_gtf):
		raise IOError("# Error: file {} does not exist".format(args.in_gtf))
	if not os.path.isfile(args.in_bam):
		raise IOError("# Error: file {} does not exist".format(args.in_bam))
	
	with open(args.out_counts,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))	

	gtf = HTSeq.GFF_Reader(str(args.in_gtf))
	bam = HTSeq.BAM_Reader(str(args.in_bam))
	
	# make a dictionary of genes with intervals that include the 5' and 3' UTRs
	print("Making Databases for {} using {}".format(args.in_bam,args.in_gtf))
	print("Exporting to {}".format(args.out_counts))
	gene_ivs = {}
	for feature in gtf:
		if feature.type == "ncRNA_gene":
			if 'Name' in feature.attr: # names seem slightly preferable to gene_ids if they exist, but all ncRNAs have gene_ids
				name = feature.attr['Name']
			else:
				name = feature.attr['gene_id']
		else:
			continue
		gene_ivs[name] = feature

	# make a feature array with the intervals of all the genes and various dictionaries
	gene_array = HTSeq.GenomicArrayOfSets("auto",stranded=True)
	counts = {}
	for gene in gene_ivs:
		counts[gene] = 0
		gene_array[gene_ivs[gene].iv] += gene

	# algorithm partially based on htseq documentation
	# https://htseq.readthedocs.io/en/master/counting.html?highlight=paired-end#handling-paired-end-reads
	# and https://htseq.readthedocs.io/en/master/tour.html
	def invert_strand(iv):
		# from HTSeq/scripts/utils.py
		iv2 = iv.copy()
		if iv2.strand == "+":
			iv2.strand = "-"
		elif iv2.strand == "-":
			iv2.strand = "+"
		else:
			raise ValueError("Illegal strand")
		return iv2

	counter = 0
	warn_sort_flag = False
	print("Counting reads:")
	# read 1 is from the opposite strand!!!
	
	test = HTSeq.GenomicInterval("IX",393844,397082)
	counter=1
	for bundle in HTSeq.pair_SAM_alignments( bam[test], bundle=True ):
		if counter % 100000 == 0:
			print(counter)
		if len(bundle) != 1:
			continue  # Skip multiple alignments
		second_almnt, first_almnt = bundle[0]  # extract pair, NEBnext puts directional strand on the second read
		if not first_almnt or not second_almnt:
			warn_sort_flag = True
			continue
		
		counter+=1
		if not first_almnt.aligned or not second_almnt.aligned:
			continue
		if args.invert_strands:
			first_almnt.iv = invert_strand(first_almnt.iv)
			second_almnt.iv = invert_strand(second_almnt.iv)
		iset = set() # find the gene that is common to both reads
		for iv, step_set in gene_array[first_almnt.iv].steps():
			if len(step_set) == 0:
				continue
			elif len(iset) == 0:
				iset = step_set.copy()
			else:
				iset.intersection_update(step_set)
		if len(iset) == 0: #move on if no genes found for the first read
			continue
		second_almnt.iv = invert_strand(second_almnt.iv) # other read is on the other strand!!
		for iv, step_set in gene_array[second_almnt.iv].steps():
			if len(step_set) > 0: 
				iset.intersection_update(step_set)
		if len(iset) != 1: # ignore ambiguous reads
			continue
		found_gene = list(iset)[0]
		counts[found_gene] += 1
	print(counter-1)
		
	
	if warn_sort_flag:
		warnings.warn("Missing read mates, please make sure the input bam is sorted by name not coordinate")
	
	print("Exporting counts!")
	# reformat distribution counts
	dist_out = []
	for gene in dist:
		df = pd.DataFrame.from_dict(dist[gene]['Count'],orient='index',columns=['Count'])
		df['Pos'] = dist[gene]['index']
		df['Region'] = dist[gene]['region']
		df = df[df.Count != 0]
		df['ORF'] = gene
		dist_out.append(df)
	pd.concat(dist_out).to_csv(dist_uncompressed,sep='\t',mode='a',index=False)
	if dist_gzip_flag:
		with open(dist_uncompressed,'rb') as src, gzip.open(dist_compressed,'wb') as dst:
			dst.writelines(src)
		os.remove(dist_uncompressed)
	
	counts_out = []
	for gene in counts:
		df = pd.DataFrame(counts[gene],index=[gene])
		#df['length_trans']=len(dist[gene]['index'])
		df['CDS_length']=dist[gene]['CDS_len']
		df['ORF'] = df.index
		counts_out.append(df)
	counts_df = pd.concat(counts_out)
	counts_df['RPK'] = (counts_df['spliced']+counts_df['other']) / (counts_df['CDS_length']/1000)
	totalRPK = counts_df['RPK'].sum()
	counts_df['TPM'] = counts_df['RPK']*1e6/totalRPK
	counts_df.drop('RPK',axis=1).to_csv(counts_uncompressed,sep='\t',index=False,mode='a')
	if counts_gzip_flag:
		with open(counts_uncompressed,'rb') as src, gzip.open(counts_compressed,'wb') as dst:
			dst.writelines(src)
		os.remove(counts_uncompressed)
