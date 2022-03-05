"""
Created on 9/27/2021
Contains helpful functions for dealing with GFF files
"""
import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import HTSeq
import itertools
import re

def gff_todict(gff_filename):
	# reads a gff file and organizes all the features by gene
	# traces parents from feature -> mRNA -> gene
	gene_dict = {}
	transcript_dict = {}  # going to store transcripts to make it easy to find their parent genes
	gff_file = HTSeq.GFF_Reader(gff_filename,end_included=True)
	print_flag = False
	for feature in gff_file:
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
				print("No name found for {}".format(feature.name))

			# now add the gene to the gene dictionary
			gene_dict[feature.name] = { 'gene_id':gene_id,'feature':feature,
			                            'transcripts':[],'5UTR':[],'3UTR':[],
			                            'exons':[],'CDSs':[] }

	# scan again for mRNAs associated with each gene
	for feature in gff_file:
		if feature.type == "mRNA":
			if not feature.name in transcript_dict:
				transcript_dict[feature.name] = feature.attr['Parent']
			try:
				gene_dict[feature.attr['Parent']]['transcripts'].append(feature)
			except:
				print("Cannot find parent {} of transcript {}".format(feature.attr['Parent'],feature.name))

	# scan for features
	def find_parent(feature):
		parent_gene = None
		try:
			parent_gene = transcript_dict[feature.attr['Parent']]
		except:
			print("Cannot find parent of {}".format(feature.name))
		return(parent_gene)

	for feature in gff_file:
		if feature.type not in ['exon','CDS','5UTR','3UTR']:
			continue
		if not find_parent(feature):
			continue
		parent_gene = find_parent(feature)
		if feature.type == 'exon':
			gene_dict[parent_gene]['exons'].append(feature)
		elif feature.type == 'CDS':
			gene_dict[parent_gene]['CDSs'].append(feature)
		elif feature.type == '5UTR':
			gene_dict[parent_gene]['5UTR'].append(feature)
		elif feature.type == '3UTR':
			gene_dict[parent_gene]['3UTR'].append(feature)

	# output a dictionary but with gene_id as keys
	out_dict = {}
	for gene in gene_dict:
		out_dict[gene_dict[gene]['gene_id']] = {
			'feature':gene_dict[gene]['feature'],
			'transcripts':gene_dict[gene]['transcripts'],
			'5UTR':gene_dict[gene]['5UTR'],
			'3UTR':gene_dict[gene]['3UTR'],
			'exons':gene_dict[gene]['exons'],
			'CDSs':gene_dict[gene]['CDSs']
		}
	return(out_dict)