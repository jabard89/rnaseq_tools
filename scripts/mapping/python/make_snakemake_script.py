"""
Created on 2/2/2023
@author: Jared Bard
Takes a csv and outputs a bash script for analyzing via snakemake
CSV must have a column of unique Sample names, and a "SRR" header that corresponds to the fastq names
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
import itertools
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate options")
    # Required arguments
    parser.add_argument(dest="sample_sheet",default=None,type=str,help="input sample sheet")
    parser.add_argument(dest="n_fastq",default=2,type=int,help="how many fastqs (1 is SE, 2 is PE)")
    parser.add_argument(dest="transcriptome",default=None,type=str,help="name of transcriptome")
    parser.add_argument(dest="output_script",default=None,type=str,help="snakemake script")
    parser.add_argument('--fastq_dir',dest="fastq_dir",default="",type=str,help="directory for fastq files")
    args = parser.parse_args()
    # args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/seq_tutorial/CHIKV/Basavappa_2022-sample.csv",
    #                           "2",
    #                           "Homo_sapiens.GRCh38.cdna.all",
    #                           "/home/jabard89/Dropbox/code_JB/repos/seq_tutorial/CHIKV/run_snakemake.sh"])
    
    for file in [args.sample_sheet]:
        if not os.path.isfile(file):
            raise IOError("# Error: file {} does not exist".format(file))

    for file in [args.output_script]:
        if os.path.isfile(file):
            raise IOError("# Error: file {} already exists, please remove it".format(file))

    if args.n_fastq not in [1,2]:
        raise IOError("n_fastq can be 1 (SE) or 2 (PE)")

    df_samples = pd.read_csv(args.sample_sheet)

    sample=df_samples.Sample.tolist()
    
    params = "--use-conda --cores all -p"

    if args.n_fastq == 1:
        fastq1 = [args.fastq_dir+i+".fastq.gz" for i in df_samples.SRR.tolist()]
        output_tuples = list(zip(sample,fastq1))
        se_string = "snakemake --config sample={} fastq1={} transcriptome={}"+" "+params
        output_list = [se_string.format(*i,args.transcriptome) for i in output_tuples]
    elif args.n_fastq == 2:
        fastq1 = [args.fastq_dir+i+"_1.fastq.gz" for i in df_samples.SRR.tolist()]
        fastq2 = [args.fastq_dir+i+"_2.fastq.gz" for i in df_samples.SRR.tolist()]
        output_tuples = list(zip(sample,fastq1,fastq2))
        se_string = "snakemake --config sample={} fastq1={} fastq2={} transcriptome={}"+" "+params
        output_list = [se_string.format(*i,args.transcriptome) for i in output_tuples]

    output = ['#!/bin/bash'] + output_list
    with open(args.output_script, 'w') as f:
        for line in output:
            f.write(f"{line}\n")

