#!/bin/bash
curl -C - -o ./Saccharomyces_cerevisiae.R64-1-1.109.gff3.gz https://ftp.ensembl.org/pub/release-109/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gff3.gz
gzip -dfc ./Saccharomyces_cerevisiae.R64-1-1.109.gff3.gz > ./Saccharomyces_cerevisiae.R64-1-1.109.gff3
curl -C - -o ./Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gzip -dfc ./Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz > ./Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

curl -C - -o gffread-0.12.7.Linux_x86_64.tar.gz http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz
samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
tar -xvzf gffread-0.12.7.Linux_x86_64.tar.gz

gffread-0.12.7.Linux_x86_64/gffread Saccharomyces_cerevisiae.R64-1-1.109.gff3 \
	-g Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
	--w-add 600 -w Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_600down.fasta

python ./trim_fasta.py Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_600down.fasta \
	Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_30down.fasta \
	-f 0 -b 570 -m 0

python ./rename_fasta_filter_mRNA.py Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_30down.fasta \
	Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_30down_renamed.fasta