#!/bin/bash
# Download the ncRNA sequences from Ensembl
curl -C - -o ./Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz
gzip -dfc ./Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz > ./Saccharomyces_cerevisiae.R64-1-1.ncrna.fa
# list of SGD ids for rRNA
#rRNA_ids = "YNCL0016C YNCL0025C YNCL0012C YNCL0021C YNCL0018W YNCL0027W YNCL0028W YNCL0029W YNCL0030W YNCL0031W YNCL0014C YNCL0023C YNCQ0002W YNCQ0006W"
rRNA_ids="RDN18-1 RDN18-2 RDN25-1 RDN25-2 RDN5-1 RDN5-2 RDN5-3 RDN5-4 RDN5-5 RDN5-6 RDN58-1 RDN58-2 Q0020 Q0158"
ensembl_rRNA_ids=""
for id in $rRNA_ids; do
  ensembl_rRNA_ids+="${id}_rRNA "
done

# use samtools to extract the rRNA sequences
for i in $ensembl_rRNA_ids ; do 
    samtools faidx ./Saccharomyces_cerevisiae.R64-1-1.ncrna.fa ${ensembl_rRNA_ids} > Saccharomyces_cerevisiae.R64-1-1.rRNA.fa
done