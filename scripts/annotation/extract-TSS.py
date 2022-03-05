# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:43:59 2020

@author: Jared Bard
Script for pulling consensus transcript start sites for all genes from yeasTSS.org
Currently, it skips the mitochondrial chromosome
Requires a gff and fasta file for the organism, in addition to the yeasTSS consensus promoter file
"""

#Plan: From the consensus clusters, pull the dominant CTSS and the TPM of the dominant CTSS
#Read in the GFF, make a dictionary of all genes
#Also make a dictionary for each strand of each chromosome, containing all of the annotated genes and their start sites
#Match each dominant CTSS to a gene
#Pull the sequence (CTSS to CDS start) from the fasta sequence and add to gene dictionary
#Export gene dictionary (ORF, 5' UTR, TPM of dominant CTSS, organism, growth_condition)
#Then, in R, read these tables in, and cluster genes into groups using OrthoDB list
import sys, os, math, random, argparse, csv
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import gffutils
import csv
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generate ")
	# Required arguments
    parser.add_argument(dest="gff_file", default=None, type=str,help="input GFF sequence")
    parser.add_argument(dest="fasta_file", default=None, type=str,help="input FASTA sequence")
    parser.add_argument(dest="TSS_file", default=None, type=str,help="input consensus cluster file from yeastTSS")
    parser.add_argument(dest="output_file", default=None, type=str,help="file to export UTRs to")
    parser.add_argument("--organism",dest='organism',default='NA',type=str,help='organism')
    parser.add_argument("--condition",dest='condition',default='NA',type=str,help='growth condition')
    parser.add_argument("--max",dest="max_UTR",default=600,type=int,help="maximum length of UTR") #GCN4 is 572
    args = parser.parse_args()
    

    #load the GFF file
    #Extract genes to dictionary, and also create a list of the start sites for each gene, organized by chromosome
    db = gffutils.create_db(args.gff_file,':memory:',force=True,keep_order=True,merge_strategy='merge',sort_attribute_values=True)

    gene_dict = {}
    chr_dict={}
    print_flag=None
    for gene in db.features_of_type('gene',order_by='start'):
        try:
            if 'gene_id' in gene.attributes: #this is right for cerevisiae
                gene_id=gene['gene_id'][0]
                if not print_flag:
                    print("Names from the 'gene_id' attribute") #print this to the output once
                    print_flag=True
            else:
                gene_id=gene['Name'][0] #this seems better for all other provided gffs in yeasTSS
                if not print_flag:
                    print("Name from the 'Name' attribute")
                    print_flag=True
        except: print("No name found for {}".format(gene['ID']))
        #first make a new chromosome in chr_dict if necessary
        if not gene.seqid in chr_dict:
            chr_dict[gene.seqid]={'+':{'genes':[],'gene_starts':[]},
                    '-':{'genes':[],'gene_starts':[]}}
        all_CDS = [CDS for CDS in db.children(gene,featuretype='CDS',order_by='start')]
        #if the gene is not a protein, it shouldn't have any CDS features? I also tried using the biotype field, but not all of the gffs had biotypes listed.
		#currently i also pull CDSs for non-proteins, but presumably these won't have strong TSS signals
        if not all_CDS: continue
        if gene.strand=='+':
            CDS_start=all_CDS[0].start
            chr_dict[gene.seqid]['+']['genes'].append(gene_id)
            chr_dict[gene.seqid]['+']['gene_starts'].append(CDS_start)
        else: #minus strand, pull the last CDS listed
            CDS_start=all_CDS[-1].stop
            chr_dict[gene.seqid]['-']['genes'].append(gene_id)
            chr_dict[gene.seqid]['-']['gene_starts'].append(CDS_start)
        gene_dict[gene_id]={'gene_id':gene_id,'strand':gene.strand,
                 'CDS_start':CDS_start,'chr':gene.seqid,'tss':{'sites':[],'tpms':[]},
                 'main_tss':None,'main_tss_tpm':None,'5UTR':None,'5UTR_length':None}
    #load the TSS dictionary from the Yeastss file   
    TSS_dict={}
    with open(args.TSS_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader,None) #skip the header
        for row in reader:
            row_dict={'chr':row[1],'strand':row[4],'ctss':int(row[5]),'tpm':float(row[7])}
            TSS_dict[row[0]]=row_dict
    
    #Load the fasta file by chromosome
    fasta_dict={}
    for record in SeqIO.parse(args.fasta_file,format='fasta'):
        fasta_dict[record.id]=record
        
    # Go through all of the TSS and match them to the closest gene
    for TSS in TSS_dict.values():
        search_chr = TSS['chr']
        search_strand = TSS['strand']
        search_start = TSS['ctss']
        if search_chr not in chr_dict:
            continue
        strand_genes = chr_dict[search_chr][search_strand]['genes']
        strand_starts = chr_dict[search_chr][search_strand]['gene_starts']
        found_start_list = []
        found_start = None
        for start in strand_starts:
            if search_strand == '+':
                if args.max_UTR >= start-search_start > 0:
                    found_start_list.append(start)
                if found_start_list:
                    found_start = min(found_start_list)
            # opposite search for the minus strand
            else:
                if args.max_UTR >= search_start-start > 0:
                    found_start_list.append(start)
                if found_start_list:
                    found_start = max(found_start_list)
        if found_start:
            found_gene = strand_genes[strand_starts.index(found_start)]
            gene_dict[found_gene]['tss']['sites'].append(search_start)
            gene_dict[found_gene]['tss']['tpms'].append(TSS['tpm'])
    
    # find the max TPM for each gene, then extract the 5' UTR for that TSS
    # This heuristic could be changed to choose the closest
    genes_out_list=[]
    for gene in gene_dict.values():
        tss_list = gene['tss']['tpms']
        if tss_list:
            genes_out_list.append(gene['gene_id'])
            gene['main_tss']=gene['tss']['sites'][tss_list.index(max(tss_list))] #gene['tss']['sites'] and gene['tss']['tpms'] have matched indices
            gene['main_tss_tpm']=max(tss_list)
            chromo=gene['chr']
            if not chromo in fasta_dict:
                print("Chromosome {} not found".format(chromo))
                continue
            tss=gene['main_tss']
            start=gene['CDS_start']
            if gene['strand']=='+':
                gene['5UTR']=str(fasta_dict[chromo][tss-1:start-1].seq)
            else:
                gene['5UTR']=str(fasta_dict[chromo][start:tss].seq.reverse_complement())
            gene['5UTR_length']=len(gene['5UTR'])
    table_header=['ORF','main_tss','main_tss_tpm','5UTR','5UTR_length',
                   'organism','growth_condition']
    table_out=[]
    for gene_id in genes_out_list:
        gene=gene_dict[gene_id]
        table_out.append([gene_id,gene['main_tss'],gene['main_tss_tpm'],
             gene['5UTR'],gene['5UTR_length'],
             args.organism,args.condition])
    with open(args.output_file,'w',newline='') as csvfile:
        csvwriter=csv.writer(csvfile)
        csvwriter.writerow(table_header)
        csvwriter.writerows(table_out)