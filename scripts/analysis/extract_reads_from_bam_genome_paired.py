#! python
# Extracts reads mapped to every gene from data that was mapped to genomes
# Jared Bard
# 03/03/2022

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq
import re
import itertools

# reads gtf output by 

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest="in_gtf",default=None,help="input gtf file")
	parser.add_argument(dest="in_bam", default=None, help="input bam file")
	parser.add_argument(dest="out_dist",default=None,help="name of the gene distribution output tsv")
	parser.add_argument(dest="out_counts",default=None,help="name of the gene count output tsv")
	# Optional arguments
	args = parser.parse_args()
#	args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/annotations/Scerevisiae.R64-1-1.104.yeastss.pelechano.gtf",
#							"/home/jabard89/Dropbox/code_JB/repos/Polysome_seq/big/output_220303/JB049/211216_STAR_JB049_Aligned.out.bam",
#							"/home/jabard89/Dropbox/code_JB/repos/Polysome_seq/big/output_220303/JB049/211216_STAR_JB049_dist_test.tsv",
#							"/home/jabard89/Dropbox/code_JB/repos/Polysome_seq/big/output_220303/JB049/211216_STAR_JB049_counts_test.tsv"])

	# Read input
	if not os.path.isfile(args.in_gtf):
		raise IOError("# Error: file {} does not exist".format(args.in_gtf))
	if not os.path.isfile(args.in_bam):
		raise IOError("# Error: file {} does not exist".format(args.in_bam))
	with open(args.out_dist,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))
	with open(args.out_counts,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))	

	gtf = HTSeq.GFF_Reader(args.in_gtf)
	bam = HTSeq.BAM_Reader(args.in_bam)
	
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

	#extract the total interval for each gene, including 3' and 5' with a 100bp cushion
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
		if tempStart >= 100: #make sure the padding doesn't go negative
			tempStart = tempStart-100			
		gene_ivs[gene] = HTSeq.GenomicInterval(tempChrom,tempStart,tempEnd+100,tempStrand)
	
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
		pad = 100
		if len(gene_dict[gene]['5UTR']) == 1:
			iv1_start = gene_dict[gene]['5UTR'][0].iv.start_d
			iv1_end = gene_dict[gene]['5UTR'][0].iv.end_d
		else:
			iv1_start = gene_ivs[gene].start_d
			iv1_end = gene_dict[gene]['gene'][0].iv.start_d
		if strand == '-':
			transcript_parts.append(range(iv1_start+pad,iv1_end,-1))
		else:
			transcript_parts.append(range(iv1_start-pad,iv1_end))

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
			transcript_parts.append(range(iv3_start,iv3_end-100,-1))
		else:
			transcript_parts.append(range(iv3_start,iv3_end+100))

		transcript_pos = [l for part in transcript_parts for l in list(part)] #turn the ranges into a continuous vector of positions
		CDS_len = sum([len(part) for part in transcript_parts[1:-1]]) #calculate the length of the CDS
		transcript_index = list(range(0-len(transcript_parts[0]),CDS_len+len(transcript_parts[-1])))
		region_vector = ['UTR5']*len(transcript_parts[0]) + ['CDS']*CDS_len + ['UTR3']*len(transcript_parts[-1])
		dist[gene] = {'index':transcript_index,'Count':dict.fromkeys(transcript_pos,0),'region':region_vector}
		
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
	print("Counting reads:")
	# read 1 is from the opposite strand!!!
	for bundle in HTSeq.pair_SAM_alignments( bam, bundle=True ):
		counter += 1
		if counter % 100000 == 0:
			print(counter)
		if len(bundle) != 1:
			continue  # Skip multiple alignments
		second_almnt, first_almnt = bundle[0]  # extract pair
		if not first_almnt.aligned or not second_almnt.aligned:
			continue
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
		# now check if any read end is in an intron
		# using this cigar loop from
		# https://htseq.readthedocs.io/en/master/counting.html?highlight=intron#cigar-operations
		# to identify any subregion of the alignment which may match to an intron
		# still possible to miss introns in the gaps between the reads, but the reads are long enough relative to the fragment sizes
		# that this should be a relatively rare occurence
		feature_set = set()
		for almnt in (first_almnt,second_almnt):
			for cigop in almnt.cigar:
				if cigop.type != "M":
					continue
				for iv, val in gene_array_dict[found_gene][cigop.ref_iv].steps():
					feature_set |= val
		if 'intron' in feature_set and len(feature_set) == 1:
			counts[found_gene]['intron'] += 1
		elif 'intron' in feature_set and 'CDS' in feature_set:
			counts[found_gene]['unspliced'] += 1
		else:
			counts[found_gene]['other'] += 1
			# assign the read to the correct position on the gene
			read_pos = first_almnt.iv.start_d
			if read_pos not in dist[found_gene]['Count']:
				if first_almnt.iv.strand == "+":
					read_pos = min(dist[found_gene]['Count'].keys())
				else:
					read_pos = max(dist[found_gene]['Count'].keys())
			dist[found_gene]['Count'][read_pos] += 1
	
	print("Exporting counts!")
	# reformat distribution counts
	dist_out = []
	for gene in dist:
		df = pd.DataFrame.from_dict(dist[gene]['Count'],orient='index',columns=['Count'])
		#df['Pos'] = df.index can include position if wanted
		df['index'] = dist[gene]['index']
		df['Region'] = dist[gene]['region']
		df.set_index('index',drop=True,inplace=True)
		df = df[df.Count != 0]
		df['ORF'] = gene
		dist_out.append(df)
	pd.concat(dist_out).to_csv(args.out_dist,sep='\t',mode='a')
	
	counts_out = []
	for gene in counts:
		df = pd.DataFrame(counts[gene],index=[gene])
		df['ORF'] = df.index
		counts_out.append(df)
	pd.concat(counts_out).to_csv(args.out_counts,sep='\t',index=False,mode='a')

