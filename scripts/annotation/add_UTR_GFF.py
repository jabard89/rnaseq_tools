"""
Created on 6/11/2021
Modified 9/26/2021
@author: Jared Bard
Script for pulling consensus transcript start sites for all genes from yeasTSS.org
and then integrating them into a GTF file (which is a subversion of a GFF file)

uses HTseq to read and write the GTF file
Requires a gff for the organism, in addition to the yeasTSS consensus promoter file

Based on hs-datasets/utility/YeasTss/extract-TSS.py

Reading in latest GFF from ensembl, also 3' UTRs from Pelechano 2013
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import HTSeq
import itertools
import re

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF sequence")
	parser.add_argument(dest="TSS_file",default=None,type=str,help="input consensus cluster file from yeastTSS")
	parser.add_argument(dest="Pelechano_file",default=None,type=str,
	                    help="input Pelechano 3' UTRs downloaded from yeastgenome.org")
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export GTF to")
	parser.add_argument("--max",dest="max_UTR",default=600,type=int,help="maximum length of UTR")  # GCN4 is 572
	args = parser.parse_args()
	# args = parser.parse_args(["src/Saccharomyces_cerevisiae.R64-1-1.104.gff3",
	#                           "src/Scerevisiae/ScerYPDconsensusClusters.txt",
	#                           "src/SGD_all_ORFs_3prime_UTRs.fsa",
	#                           "output/Scerevisiae.R64-1-1.104.yeastss.pelechano.gtf")

	# Write out parameters
	with open(args.output_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))

	# read the GTF file, grouping together all the CDS by gene
	gtf_file = HTSeq.GFF_Reader(args.gff_file,end_included=True)
	chr_dict = { }  # going to store chromosomes -> strands -> genes
	gene_dict = { }
	transcript_dict = { }  # going to store transcripts to make it easy to find their parent genes
	print_flag = None

	# first scan through the gtf_file and pull all of the genes
	for feature in gtf_file:
		if feature.type == "gene":
			try:
				if 'gene_id' in feature.attr:  # this is right for cerevisiae
					gene_id = feature.attr['gene_id']
					if not print_flag:
						print("Names from the 'gene_id' attribute")  # print this to the output once
						print_flag = True
				else:
					gene_id = feature.attr['Name']  # this seems better for all other provided gffs in yeasTSS
					if not print_flag:
						print("Name from the 'Name' attribute")
						print_flag = True
			except:
				print("No name found for {}".format(gene['ID']))

			# now add the gene to the gene dictionary
			gene_dict[feature.name] = { 'gene_id':gene_id,'feature':feature,
			                            'start_codon':None,'tss':{ 'sites':[],'tpms':[] },
			                            'main_tss':None,'main_tss_tpm':None,
			                            '5UTR':None,'5UTR_length':None,
			                            'transcripts':[],'exons':[],'CDSs':[]
			                            }

	# scan for mRNAs associated with each gene
	for feature in gtf_file:
		if feature.type == "mRNA":
			if not feature.name in transcript_dict:
				transcript_dict[feature.name] = feature.attr['Parent']
			try:
				gene_dict[feature.attr['Parent']]['transcripts'].append(feature)
			except:
				print("Cannot find parent {} of transcript {}".format(feature.attr['Parent'],feature.name))
	# scan for exons
	for feature in gtf_file:
		if feature.type == 'exon':
			try:
				parent_gene = transcript_dict[feature.attr['Parent']]
			except:
				print("Cannot find parent {} of exon {}".format(feature.attr['Parent'],feature.name))
			gene_dict[parent_gene]['exons'].append(feature)
	# scan for CDS
	for feature in gtf_file:
		if feature.type == 'CDS':
			try:
				parent_gene = transcript_dict[feature.attr['Parent']]
			except:
				print("Cannot find parent {} of CDS {}".format(feature.attr['Parent'],feature.name))
			gene_dict[parent_gene]['CDSs'].append(feature)

	# now find the first nucleuotide of the first CDS
	for gene in gene_dict:
		start_codon = None
		strand = gene_dict[gene]['feature'].iv.strand
		for CDS in gene_dict[gene]['CDSs']:
			CDS_start = CDS.iv.start_d + 1  # HTseq is 0-based
			if not start_codon:
				start_codon = CDS_start
			if strand == "+":  # scan for the first start codon
				if CDS_start < start_codon:
					start_codon = CDS_start
			if strand == "-":
				if CDS_start > start_codon:
					start_codon = CDS_start
		gene_dict[gene]['start_codon'] = start_codon

	# create a dictionary with start_codons
	start_dict = { }
	for gene in gene_dict:
		chrom = gene_dict[gene]['feature'].iv.chrom
		strand = gene_dict[gene]['feature'].iv.strand
		if not chrom in start_dict:
			start_dict[chrom] = { '-':{ },'+':{ } }
		if gene_dict[gene]['start_codon']:
			start_dict[chrom][strand][gene_dict[gene]['start_codon']] = gene

	TSS_dict = { }
	with open(args.TSS_file,'r') as f:
		reader = csv.reader(f,delimiter='\t')
		next(reader,None)  # skip the header
		for row in reader:
			row_dict = { 'chr':row[1],'strand':row[4],'ctss':int(row[5]),'tpm':float(row[7]) }
			TSS_dict[row[0]] = row_dict

	# Go through all of the TSS and match them to the closest gene
	found_genes_counter = 0
	for TSS in TSS_dict.values():
		search_chr = TSS['chr'][3:]
		search_strand = TSS['strand']
		search_start = TSS['ctss']
		if search_chr not in start_dict:
			# print('Missing chromomosome '+search_chr)
			continue
		strand_starts = start_dict[search_chr][search_strand].keys()
		found_start_list = []
		found_start = None
		for start in strand_starts:
			if search_strand == '+':
				if args.max_UTR >= start - search_start > 0:
					found_start_list.append(start)
				if found_start_list:
					found_start = min(found_start_list)
			# opposite search for the minus strand
			else:
				if args.max_UTR >= search_start - start > 0:
					found_start_list.append(start)
				if found_start_list:
					found_start = max(found_start_list)
		if found_start:
			found_genes_counter += 1
			found_gene = start_dict[search_chr][search_strand][found_start]
			gene_dict[found_gene]['tss']['sites'].append(search_start)
			gene_dict[found_gene]['tss']['tpms'].append(TSS['tpm'])
	print("Found {} genes matched with TSS".format(found_genes_counter))

	# find the max TPM for each gene, then extract the 5' UTR for that TSS
	# This heuristic could be changed to choose the closest
	genes_out_list = []
	for gene in gene_dict.values():
		tss_list = gene['tss']['tpms']
		if tss_list:
			genes_out_list.append(gene['gene_id'])
			gene['main_tss'] = gene['tss']['sites'][
				tss_list.index(max(tss_list))]  # gene['tss']['sites'] and gene['tss']['tpms'] have matched indices
			gene['main_tss_tpm'] = max(tss_list)
			tss = gene['main_tss']
			start = gene['start_codon']
			if gene['feature'].iv.strand == '+':
				gene['5UTR_start'] = tss
				gene['5UTR_end'] = start - 1
				gene['5UTR_length'] = abs(start - tss) - 1
			else:
				gene['5UTR_start'] = start + 1
				gene['5UTR_end'] = tss
				gene['5UTR_length'] = abs(start - tss)

	# now load the 3' UTRs
	with open(args.Pelechano_file,'r') as f:
		Pelechano = [line for line in f if line[0] == '>']
	Pelechano_dict = { }
	for line in Pelechano:
		ORF_search = re.search('(?<=_3950_).*(?=_id)',line)
		if ORF_search:
			gene_id = ORF_search.group(0)
		else:
			next
		if not gene_id in Pelechano_dict:
			Pelechano_dict[gene_id] = { 'gene_id':gene_id,'chr':[],'start':[],'stop':[],'length':[] }

		range_search = re.search('(?<=range=chr).*(?= 5)',line)
		if range_search:
			range1 = range_search.group(0)
		else:
			print("Could not find range for {}".format(gene_id))
			break

		strand_search = re.search('(?<=strand=).',line)
		if strand_search:
			strand = strand_search.group(0)

		chr_search = re.search('.*(?=:)',range1)
		if chr_search:
			chr = chr_search.group(0)
			Pelechano_dict[gene_id]['chr'].append(chr)
		else:
			print("Could not find chr for {}".format(gene_id))
			break

		start_search = re.search('(?<=:).*(?=-)',range1)
		if start_search:
			if strand == "+":
				start = int(start_search.group(
					0)) + 1  # database includes the last nt of the stop codon, correcting for that
			else:
				start = int(start_search.group(0))
			Pelechano_dict[gene_id]['start'].append(start)
		else:
			print("Could not find start for {}".format(gene_id))
			break

		stop_search = re.search('(?<=-).*',range1)
		if stop_search:
			if strand == "-":
				stop = int(stop_search.group(
					0)) - 1  # database includes the last nt of the stop codon, correcting for that
			else:
				stop = int(stop_search.group(0))
			Pelechano_dict[gene_id]['stop'].append(stop)
		else:
			print("Could not find stop for {}".format(gene_id))
			break
		Pelechano_dict[gene_id]['length'].append(abs(start - stop) + 1)
	for gene_id in Pelechano_dict.keys():
		max_length = 0
		for i,length in enumerate(Pelechano_dict[gene_id]['length']):
			if length > max_length:
				Pelechano_dict[gene_id]['3UTR_start'] = Pelechano_dict[gene_id]['start'][i]
				Pelechano_dict[gene_id]['3UTR_stop'] = Pelechano_dict[gene_id]['stop'][i]
				Pelechano_dict[gene_id]['3UTR_length'] = Pelechano_dict[gene_id]['length'][i]
				max_length = length

	# add new features to the GTF file!!
	# going to iteratre through the GTF file. For every CDS, look to see if I can find a UTR5 or UTR3 from the relevent dictionaries to add
	# meanwhile print out a new list of the gtf for rewriting
	# remember: HTSeq.GenomicInterval is 0-based, but TSS etc are stored as 1-based.
	# remember: HTSeq.GenomicInterval doesn't include the end (python style)
	gtf_file = HTSeq.GFF_Reader(args.gff_file,end_included=True)
	gtf_out = []
	for i,feature in enumerate(gtf_file):
		feature_3UTR = ""
		feature_5UTR = ""
		if feature.type == 'gene':
			gene_id = feature.attr['gene_id']
			strand = feature.iv.strand
			chrom = feature.iv.chrom
			if feature.name in gene_dict:
				g = gene_dict[feature.name]
				if g['main_tss']:
					interval = HTSeq.GenomicInterval(chrom,g['5UTR_start'] - 1,g['5UTR_end'],strand)
					feature_5UTR = HTSeq.GenomicFeature(gene_id,'5UTR',interval)
					feature_5UTR.source = 'YeasTSS'
					if g['transcripts']:
						feature_5UTR.attr = { 'ID':"5UTR:" + gene_id,'parent':g['transcripts'][0].name,
						                      'protein_id':gene_id
						                      }
					else:
						feature_5UTR.attr = { 'ID':"5UTR:" + gene_id,'parent':"transcript:" + gene_id,
						                      'protein_id':gene_id
						                      }
			if gene_id in Pelechano_dict:
				g = Pelechano_dict[gene_id]
				interval = HTSeq.GenomicInterval(chrom,g['3UTR_start'] - 1,g['3UTR_stop'],strand)
				feature_3UTR = HTSeq.GenomicFeature(gene_id,'3UTR',interval)
				feature_3UTR.source = 'Pelechano_2013'
				gene_name = 'gene:' + gene_id
				if gene_dict[gene_name]['transcripts']:
					feature_3UTR.attr = { 'ID':"3UTR:" + gene_id,'parent':gene_dict[gene_name]['transcripts'][0].name,
					                      'protein_id':gene_id
					                      }
				else:
					feature_3UTR.attr = { 'ID':"3UTR:" + gene_id,'parent':"transcript:" + gene_id,'protein_id':gene_id }
			feature.attr['description'] = ""
			gtf_out.append(feature.get_gff_line())
			if feature_5UTR: gtf_out.append(feature_5UTR.get_gff_line())
			if feature_3UTR: gtf_out.append(feature_3UTR.get_gff_line())
		else:
			feature.attr['description'] = ""
			gtf_out.append("###\n")
			gtf_out.append(feature.get_gff_line())

	# write the file!
	with open(args.output_file,'a') as f:
		f.writelines(gtf_out)

