#! python
# Extracts reads mapped to every gene from data that was mapped to genomes
# Jared Bard
# 03/03/2022
# Modified 04/05/2022 to pile up reads rather than count starts and to export TPMS from the non-intron reads
# can output either tsv or tsv.gz files

import sys, os, math, random, argparse, csv, warnings, gzip, pathlib, re, itertools
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument('-i',"--invert_strands",action="store_true",help="invert the strands when looking for matches")
	parser.add_argument(dest="in_gtf",type=pathlib.Path,default=None,help="input gtf file")
	parser.add_argument(dest="in_bam", type=pathlib.Path,default=None,help="input bam file")
	parser.add_argument(dest="out_dist",type=pathlib.Path,default=None,help="name of the gene distribution output tsv")
	parser.add_argument(dest="out_counts",type=pathlib.Path,default=None,help="name of the gene count output tsv")
	# Optional arguments
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/annotations/Scerevisiae.R64-1-1.104.yeastss.pelechano.gtf",
	# 						"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/220301/star/HG49/HG49_YFL039C.bam",
	# 						"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/220301/star/HG49/HG49_dist_test.tsv.gz",
	# 						"/home/jabard89/Dropbox/Data_JB/RNA-seq/JB/220301/star/HG49/HG49_counts_test.tsv"])

	# Read input
	if not os.path.isfile(args.in_gtf):
		raise IOError("# Error: file {} does not exist".format(args.in_gtf))
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

	gtf = HTSeq.GFF_Reader(str(args.in_gtf))
	bam = HTSeq.BAM_Reader(str(args.in_bam))
	
	# make a dictionary of genes with intervals that include the 5' and 3' UTRs
	print("Making Databases for {} using {}".format(args.in_bam,args.in_gtf))
	print("Exporting to {} and {}".format(args.out_dist,args.out_counts))
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

	#extract the total interval for each gene, including 3' and 5' with a cushion
	pad = 50
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
		if tempStart >= pad: #make sure the padding doesn't go negative
			tempStart = tempStart-pad			
		gene_ivs[gene] = HTSeq.GenomicInterval(tempChrom,tempStart,tempEnd+pad,tempStrand)
	
	# make a feature array with the intervals of all the genes and various dictionaries
	gene_array = HTSeq.GenomicArrayOfSets("auto",stranded=True)
	gene_array_dict = {} # a GenomicArray for every individual gene
	counts = {}
	dist = {}
	for gene in gene_ivs:
		counts[gene] = {'intron':0,'unspliced':0,'other':0}
		gene_array[gene_ivs[gene]] += gene
		gene_array_dict[gene] = HTSeq.GenomicArrayOfSets("auto",stranded=True)
		for feature in gene_dict[gene]:
			if feature in ['CDS','5UTR','3UTR']:
				for item in gene_dict[gene][feature]:
					gene_array_dict[gene][item.iv] += feature

		# now define introns as regions within the original gene boundaries but not containing CDS
		CDS_interval = gene_dict[gene]['gene'][0].iv
		intron_ivs = []
		for iv,step_set in gene_array_dict[gene][CDS_interval].steps():
			if len(step_set) == 0:
				intron_ivs.append(iv)
				gene_dict[gene]['intron'].append(HTSeq.GenomicFeature('intron:'+gene,'intron',iv))
		for iv in intron_ivs:
			gene_array_dict[gene][iv] += 'intron'
		
		#make the distribution dictionary by assembling all UTRs and CDS into one mapping
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

		transcript_pos = [l for part in transcript_parts for l in list(part)] #turn the ranges into a continuous vector of positions
		CDS_len = sum([len(part) for part in transcript_parts[1:-1]]) #calculate the length of the CDS
		transcript_index = list(range(0-len(transcript_parts[0]),CDS_len+len(transcript_parts[-1])))
		region_vector = ['UTR5']*len(transcript_parts[0]) + ['CDS']*CDS_len + ['UTR3']*len(transcript_parts[-1])
		dist[gene] = {'index':transcript_index,'Count':dict.fromkeys(transcript_pos,0),'region':region_vector,'CDS_len':CDS_len}
		
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
	for read in bam:
		counter += 1
		if counter % 100000 == 0:
			print(counter)
		if not read.aligned:
			continue
		if args.invert_strands:
			read.iv = invert_strand(read.iv)
		iset = set() # find the gene that is common to both reads
		for iv, step_set in gene_array[read.iv].steps():
			if len(step_set) == 0:
				continue
			elif len(iset) == 0:
				iset = step_set.copy()
			else:
				iset.intersection_update(step_set)
		if len(iset) != 1: # ignore ambiguous reads
			continue
		found_gene = list(iset)[0]
		# now check if any read end is in an intron
		# using this cigar loop from
		# https://htseq.readthedocs.io/en/master/counting.html?highlight=intron#cigar-operations
		# to identify any subregion of the alignment which may match to an intron
		# still possible to miss introns in the gaps between the reads, but the reads are long enough relative to the fragment sizes
		# that this should be a relatively rare occurence
		feature_set = set()
		for cigop in read.cigar:
			if cigop.type != "M":
				continue
			for iv, val in gene_array_dict[found_gene][cigop.ref_iv].steps():
				feature_set |= val
		if 'intron' in feature_set and len(feature_set) == 1:
			counts[found_gene]['intron'] += 1
		elif 'intron' in feature_set and 'CDS' in feature_set:
			counts[found_gene]['unspliced'] += 1
		else:
			if 'CDS' in feature_set:
				counts[found_gene]['other'] += 1
			# assign the read to the correct position on the gene
			# construct a new interval that spans the whole read (from first base of alnmt1 to last base of almnt 2)
			# add 1 to every position within that interval
			range_start = read.iv.start_d
			range_end = read.iv.end_d
			# squish down reads to fit the confines of the gene vector
			if range_start not in dist[found_gene]['Count']:
				if read.iv.strand == "+":
					range_start = min(dist[found_gene]['Count'].keys())
				else:
					range_start = max(dist[found_gene]['Count'].keys())
			if range_end not in dist[found_gene]['Count']:
				if read.iv.strand == "+":
					range_end = max(dist[found_gene]['Count'].keys())
				else:
					range_end = min(dist[found_gene]['Count'].keys())
			if range_start < range_end: # need to make sure the range goes small to large
				count_range = range(range_start,range_end+1) #+1 to account for python ranges not including last number
			else:
				count_range = range(range_end,range_start+1)
			for pos in count_range:
				if pos in dist[found_gene]['Count']:
					dist[found_gene]['Count'][pos] += 1
	
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
	counts_df['RPK'] = counts_df['other'] / (counts_df['CDS_length']/1000)
	totalRPK = counts_df['RPK'].sum()
	counts_df['TPM'] = counts_df['RPK']*1e6/totalRPK
	counts_df.drop('RPK',axis=1).to_csv(counts_uncompressed,sep='\t',index=False,mode='a')
	if counts_gzip_flag:
		with open(counts_uncompressed,'rb') as src, gzip.open(counts_compressed,'wb') as dst:
			dst.writelines(src)
		os.remove(counts_uncompressed)