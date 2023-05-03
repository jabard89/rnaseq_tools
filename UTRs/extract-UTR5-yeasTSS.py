# -*- coding: utf-8 -*-
"""
Created 230424
@author: Jared Bard
Uses raw CTSS data from yeastTSS.org to generate a table of UTRs by selecting the median
CTSS for each gene and then pulling the sequence from the fasta file
Modified from extract-TSS.py
Plan:
make a dictionary of genes from a gff file, and the 700nt upstream of that gene
remove nts that overlap with other genes
extract the CTSS values from the yeastTSS file (merge two bioreps)
find the median CTSS for each gene
pull the sequence from the fasta file
"""

import sys, os, math, random, argparse, csv,itertools
from datetime import datetime
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import csv
import HTSeq
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generate ")
    # Required arguments
    parser.add_argument(dest="gff_file", default=None, type=str,help="input GFF sequence")
    parser.add_argument(dest="fasta_file", default=None, type=str,help="input FASTA sequence")
    parser.add_argument(dest="CTSS_file1", default=None, type=str,help="input CTSS file1 from yeastTSS")
    parser.add_argument(dest="CTSS_file2", default=None, type=str,help="input CTSS file2 from yeastTSS")
    parser.add_argument(dest="output_file", default=None, type=str,help="file to export UTRs to")
    parser.add_argument("--organism",dest='organism',default='NA',type=str,help='organism')
    parser.add_argument("--condition",dest='condition',default='NA',type=str,help='growth condition')
    parser.add_argument("--max",dest="max_UTR",default=600,type=int,help="maximum length of UTR") #GCN4 is 572
    args = parser.parse_args()
    # src = "/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools"
    # args = parser.parse_args([src+"/UTRs/Saccharomyces_cerevisiae.R64-1-1.109.gff3",
    #                           src+"/UTRs/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
    #                           src+"/UTRs/yeasTSS/ScerYPD.1.ctss",
    #                           src+"/UTRs/yeasTSS/ScerYPD.2.ctss",
    #                           src+"/UTRs/yeastTSS_UTRs.csv",
    #                           "--organism=Scerevisiae","--condition=YPD","--max=700"])

    #load the GFF file
    #Extract genes to dictionary, and also create a list of the start sites for each gene, organized by chromosome
    with open(args.output_file,'w') as f:
        f.write("# Run started {}\n".format(datetime.now()))
        f.write("# Command: {}\n".format(' '.join(sys.argv)))
        f.write("# Parameters:\n")
        optdict = vars(args)
        for (k,v) in optdict.items():
            f.write("#\t{k}: {v}\n".format(k=k,v=v))
     
    
    # make a genomic array of sets to search against
    gff_HTSeq = HTSeq.GFF_Reader(args.gff_file,end_included=True)
    genes_HTSeq = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    gene_dict = {}
    gene_attr = None
    for feature in gff_HTSeq:
        if feature.type == 'gene':
            try:
                if 'gene_id' in feature.attr: #this is right for cerevisiae
                    gene_id=feature.attr['gene_id']
                    if not gene_attr:
                        gene_attr='gene_id'
                else:
                    gene_id=feature.attr['Name'] #this seems better for all other provided gffs in yeasTSS
                    if not gene_attr:
                        gene_attr='Name'
                genes_HTSeq[feature.iv] += gene_id
                gene_dict[gene_id] = {'feature':feature,'ctss':[],'TSS':None}
            except:
                print("No name found for {}".format(feature.name))
    print("Name from the {} attribute".format(gene_attr))
    
    #Load the fasta file by chromosome
    fasta_dict={}
    for record in SeqIO.parse(args.fasta_file,format='fasta'):
        fasta_dict[record.id]=record

    CTSS_1 = pd.read_csv(args.CTSS_file1,sep='\t',header=None,names=['chr','pos','strand','ctss'])
    CTSS_2 = pd.read_csv(args.CTSS_file2,sep='\t',header=None,names=['chr','pos','strand','ctss'])
    CTSS = pd.concat([CTSS_1,CTSS_2],axis=0,ignore_index=True)

    # relabel the chromosomes to match the fasta file
    if "chr" not in next(iter(fasta_dict)): 
        CTSS['chr'] = CTSS['chr'].str.replace('chr','')
    CTSS['chr'] = CTSS['chr'].str.replace('M','Mito')
    
    # extract the median CTSS for each gene
    # HTSeq is 0-based, so add 1 to the start position
    def create_temp_pos(chrom,pos,strand):
        # check if position is less than 1, if so, set to 1
        if pos <= 0:
            pos = 1
        return(HTSeq.GenomicPosition(chrom,pos-1,strand=strand)) # fixes indexing
    for gene in gene_dict.keys():
        gene_iv = gene_dict[gene]['feature'].iv
        if gene_iv.strand == "+": 
            search_end = gene_iv.start_d+1-1
            search_start = search_end-args.max_UTR+1 # +1 because closed interval
            for pos in range(search_end,search_start-1,-1):
                temp_pos = create_temp_pos(chrom=gene_iv.chrom,pos=pos,strand=gene_iv.strand)
                if len(genes_HTSeq[temp_pos]) > 0: # found an overlapping gene!
                    search_start = pos+1
                    break
        elif gene_iv.strand == "-":
            search_start = gene_iv.start_d+1+1
            search_end = search_start+args.max_UTR-1 # -1 because closed interval
            for pos in range(search_start,search_end+1,1):
                temp_pos = create_temp_pos(chrom=gene_iv.chrom,pos=pos,strand=gene_iv.strand)
                if len(genes_HTSeq[temp_pos]) > 0: # found an overlapping gene!
                    search_end = pos-1
                    break
        
        if search_start>search_end:
            gene_dict[gene]['ctss_med'] = None
            continue       
        
        # find first overlapping gene
        CTSS_array = CTSS[(CTSS['chr']==gene_iv.chrom) &
                          (CTSS['strand']==gene_iv.strand) &
                          CTSS['pos'].between(search_start,search_end,inclusive='both')]
        
        # if gene_iv.strand == "+":
        #     CTSS_array = CTSS_array.sort_values(by=['pos'],ascending=False)
        # elif gene_iv.strand == "-":
        #     CTSS_array = CTSS_array.sort_values(by=['pos'],ascending=True)
        # end_i = len(CTSS_array)
        # for i,pos in enumerate(CTSS_array['pos']):
        #     temp_pos = HTSeq.GenomicPosition(gene_iv.chrom,pos-1,strand=gene_iv.strand)
        #     if len(genes_HTSeq[temp_pos]) > 0: # found an overlapping gene!
        #         end_i = i
        # CTSS_array = CTSS_array.iloc[:end_i,:] # remove all CTSS that overlap with other genes
            
                    
        # drop_list = []
        # for index,row in CTSS_array.iterrows():
        #     temp_pos = HTSeq.GenomicPosition(gene_iv.chrom,row["pos"]-1,strand=gene_iv.strand)
        #     if len(genes_HTSeq[temp_pos]) > 0:
        #         drop_list.append(index)
        # CTSS_array = CTSS_array.drop(drop_list)

        if len(CTSS_array) == 0:
            gene_dict[gene]['ctss_med'] = None
            continue
        
        CTSS_pos = CTSS_array['pos'].to_numpy()
        CTSS_val = CTSS_array['ctss'].to_numpy()
        if gene_iv.strand == "+":
            CTSS_med = int(np.floor(np.median(np.repeat(CTSS_pos,CTSS_val))))
        elif gene_iv.strand == "-":
            CTSS_med = int(np.ceil(np.median(np.repeat(CTSS_pos,CTSS_val))))
        
        gene_dict[gene]['ctss_med'] = CTSS_med
        print(gene,CTSS_med)
    
    genes_out_list=[]
    for gene in gene_dict.keys():
        if gene_dict[gene]['ctss_med']:
            gene_iv = gene_dict[gene]['feature'].iv
            tss=gene_dict[gene]['ctss_med']
            start=gene_iv.start_d+1
            if gene_iv.strand=='+':
                gene_dict[gene]['5UTR']=str(fasta_dict[gene_iv.chrom][tss-1:start-1].seq)
            elif gene_iv.strand=='-':
                gene_dict[gene]['5UTR']=str(fasta_dict[gene_iv.chrom][start:tss].seq.reverse_complement())
            gene_dict[gene]['5UTR_length']=len(gene_dict[gene]['5UTR'])
            genes_out_list.append(gene)
    table_header=['ORF','med_tss','UTR5','UTR5_length',
                   'organism','growth_condition']
    table_out=[]
    for gene_id in genes_out_list:
        gene=gene_dict[gene_id]
        table_out.append([gene_id,gene['ctss_med'],
             gene['5UTR'],gene['5UTR_length'],
             args.organism,args.condition])
    with open(args.output_file,'a',newline='') as csvfile:
        csvwriter=csv.writer(csvfile)
        csvwriter.writerow(table_header)
        csvwriter.writerows(table_out)