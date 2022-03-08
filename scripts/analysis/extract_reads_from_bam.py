#! python

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
import HTSeq
import re

# reads gtf assuming riboviz format. may need to mess with bed reading if not in gtf format

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generic script template")
    # Required arguments
    parser.add_argument(dest="in_gtf",default=None,help="input gtf file")
    parser.add_argument(dest="in_bam", default=None, help="input bam file")
    parser.add_argument(dest="out_file",default=None,help="name of the output tsv.gz")
    # Optional arguments
    options = parser.parse_args()

    # Read input
    if not os.path.isfile(options.in_gtf):
        raise IOError("# Error: file {} does not exist".format(options.in_gtf))
    if not os.path.isfile(options.in_bam):
        raise IOError("# Error: file {} does not exist".format(options.in_bam))
    with open(options.out_file,'w') as f:
        f.write("# Run started {}\n".format(datetime.now()))
        f.write("# Command: {}\n".format(' '.join(sys.argv)))
        f.write("# Parameters:\n")
        optdict = vars(options)
        for (k,v) in optdict.items():
            f.write("#\t{k}: {v}\n".format(k=k,v=v))
            
    gtf = HTSeq.GFF_Reader(options.in_gtf)
    bam = HTSeq.BAM_Reader(options.in_bam)

    #process bam assuming that gtf is in riboviz format
    counts = {}
    for feature in gtf:
        if feature.type == "UTR3":
            if feature.name not in counts:
                counts[feature.name] = np.zeros(feature.iv.end + 1,dtype=int)
    read_count = 0
    for alngt in bam:
        read_count += 1
        if alngt.aligned:
            counts[alngt.iv.chrom][alngt.iv.start] += 1
        if (read_count%100000 == 0):
            print(read_count)
    data = []
    for ORF in counts:
        df = pd.DataFrame({"ORF":ORF,"Count":counts[ORF],"Pos":range(-250,len(counts[ORF]) - 250)})
        df = df[df.Count != 0]
        data.append(df)
    pd.concat(data).to_csv(options.out_file,sep='\t',mode='a',index=False)

