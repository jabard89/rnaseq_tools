"""
Created on 2/2/2023
@author: Jared Bard
Takes a csv and outputs a bash scripts for downloading the fastqs
CSV must have fastq urls, and a "SRR" header that corresponds to the fastq names
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
    parser.add_argument(dest="output_script",default=None,type=str,help="output script")
    args = parser.parse_args()
    # args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/seq_tutorial/CHIKV/Basavappa_2022-sample.csv",
    #                           "2",
    #                           "/home/jabard89/Dropbox/code_JB/repos/seq_tutorial/CHIKV/fastq/download_fastq.sh"])
    
    for file in [args.sample_sheet]:
        if not os.path.isfile(file):
            raise IOError("# Error: file {} does not exist".format(file))

    for file in [args.output_script]:
        if os.path.isfile(file):
            raise IOError("# Error: file {} already exists, please remove it".format(file))

    if args.n_fastq not in [1,2]:
        raise IOError("n_fastq can be 1 (SE) or 2 (PE)")

    df_samples = pd.read_csv(args.sample_sheet)

    download1 = ["wget -N "+i for i in df_samples.fastq1.tolist()]

    download2 = []

    if args.n_fastq == 2:
        download2 = ["wget -N "+i for i in df_samples.fastq2.tolist()]
    
    fastq_output = ['#!/bin/bash'] + download1 + download2
    with open(args.output_script, 'w') as f:
        for line in fastq_output:
            f.write(f"{line}\n")

