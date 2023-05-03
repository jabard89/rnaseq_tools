"""
Created on 04/24/2023
@author: Jared Bard
Based on gffutils_add_geneid.py
Takes a TSV file with UTR3_length and UTR5_length, and adds these as attributes to the GFF3 file
Also adds gene_id to each feature if it is not already present
Input file format:
ORF\tUTR5_length\tUTR3_length
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import gffutils

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF3 file")
	parser.add_argument(dest="utr_file",default=None,type=str,help="input utr file")
	parser.add_argument("-src1","--utr5_src",default="none",type=str,help="source of UTR5_length")
	parser.add_argument("-src2","--utr3_src",default="none",type=str,help="source of UTR3_length")
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export gff3 to")
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/"
	# 		   					"UTRs/Saccharomyces_cerevisiae.R64-1-1.109.gff3",  
	# 							"/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/"
	# 							"UTRs/Saccharomyces_cerevisiae_YPD_combined_yeasTSS_Pelechano2013_UTRs.tsv",                        
	#                           	"/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/"
	# 		   					"UTRs/Saccharomyces_cerevisiae.R64-1-1.109_yeasTSS_Pelechano.gff3",
	# 							"--utr5_src=Pelechano2013","--utr3_src=Pelechano2013"])
	
	for file in [args.gff_file,args.utr_file]:
		if not os.path.isfile(file):
			raise IOError("# Error: file {} does not exist".format(file))	
	
	for file in [args.output_file]:
		if os.path.isfile(file):
			raise IOError("# Error: file {} already exists, please remove it".format(file))
	
	# load UTR file
	df_utrs = pd.read_csv(args.utr_file,sep='\t',comment="#",
		       dtype={'ORF':"str",'UTR5_length':"Int64",'UTR3_length':"Int64"})
	# Write out parameters
	with open(args.output_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))
		
	db = gffutils.create_db(args.gff_file, ":memory:",
			 id_spec={"gene":"ID","mRNA":"ID","exon":"Name"})
	
	with open(args.output_file, 'a') as fout:
		for d in db.directives:
			fout.write('## {0}\n'.format(d))

		for feature in db.all_features():
			if feature.featuretype in ["exon","CDS","mRNA"]:
				feature['ID'] = feature.id
				gene_parent = list(db.parents(feature,featuretype='gene'))
				if gene_parent:
					if len(gene_parent) > 1:
						print("feature {} has multiple gene parents".format(feature.id))
					else:
						gene = gene_parent[0]['gene_id'][0]
						feature['gene_id'] = gene
				fout.write(str(feature) + '\n')

				if feature.featuretype == 'mRNA':
					parent_id = feature.id
					utr_row = df_utrs[df_utrs['ORF'] == gene].squeeze()
					utr5 = None
					utr3 = None
					if (len(utr_row)!=0) & (feature.strand == "+"):
						if not pd.isna(utr_row['UTR5_length']):
							utr5 = [feature.chrom,args.utr5_src,"five_prime_UTR",
										feature.start-utr_row['UTR5_length'],
										feature.start-1,
										".",feature.strand,".",
										"ID=five_prime_utr:{};Parent={};gene_id={}".
										format(gene,parent_id,gene)]
						if not pd.isna(utr_row['UTR3_length']):
							utr3 = [feature.chrom,args.utr3_src,"three_prime_UTR",
										feature.stop+1,
										feature.stop+utr_row['UTR3_length'],
										".",feature.strand,".",
										"ID=three_prime_utr:{};Parent={};gene_id={}".
										format(gene,parent_id,gene)]
					elif (len(utr_row!=0)) & (feature.strand == "-"):
						if not pd.isna(utr_row['UTR5_length']):
							utr5 = [feature.chrom,"yeasTSS","five_prime_UTR",
										feature.stop+1,
										feature.stop+utr_row['UTR5_length'],
										".",feature.strand,".",
										"ID=five_prime_utr:{};Parent={};gene_id={}".
										format(gene,parent_id,gene)]
						if not pd.isna(utr_row['UTR3_length']):
							utr3 = [feature.chrom,"Pelechano_2013","three_prime_UTR",
										feature.start-utr_row['UTR3_length'],
										feature.start-1,
										".",feature.strand,".",
										"ID=three_prime_utr:{};Parent={};gene_id={}".
										format(gene,parent_id,gene)]
					if utr5:
						utr5_str = '\t'.join([str(x) for x in utr5])
						fout.write(str(gffutils.feature.feature_from_line(utr5_str)) + '\n')
					if utr3:
						utr3_str = '\t'.join([str(x) for x in utr3])
						str(gffutils.feature.feature_from_line(utr3_str))
						fout.write(str(gffutils.feature.feature_from_line(utr3_str)) + '\n')
			else:
				fout.write(str(feature) + '\n')
			
		

	