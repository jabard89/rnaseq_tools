#!/bin/bash

POMBE_FA="Schizosaccharomyces_pombe_all_chromosomes.fa"
POMBE_GFF3="Schizosaccharomyces_pombe_all_chromosomes.gff3"
wget -nc -O "../annotation/${POMBE_FA}.gz" "https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
gzip -dfc ../annotation/${POMBE_FA}.gz > ../annotation/${POMBE_FA}

wget -nc -O "../annotation/${POMBE_GFF3}.gz" "https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"
gzip -dfc ../annotation/${POMBE_GFF3}.gz > ../annotation/${POMBE_GFF3}

sed -i 's/>/>pombe/g' ../annotation/${POMBE_FA}
sed -i 's/^I/pombeI/g' ../annotation/${POMBE_GFF3}
sed -i 's/^II/pombeII/g' ../annotation/${POMBE_GFF3}
sed -i 's/^III/pombeIII/g' ../annotation/${POMBE_GFF3}
sed -i 's/^MT/pombeMT/g' ../annotation/${POMBE_GFF3}
sed -i 's/^MTR/pombeMTR/g' ../annotation/${POMBE_GFF3}
sed -i 's/^AB325691/pombeAB325691/g' ../annotation/${POMBE_GFF3}

wget -nc -O ../annotation/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz http://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gzip -dfc ../annotation/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz > ../annotation/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

wget -nc -O ../annotation/sgd_20210422_orf_coding_all.fasta.gz http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz
gzip -dfc ../annotation/sgd_20210422_orf_coding_all.fasta.gz > ../annotation/sgd_20210422_orf_coding_all.fasta

wget -nc -O ../annotation/sgd_20191025_rna_coding.fasta.gz http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz
gzip -dfc ../annotation/sgd_20191025_rna_coding.fasta.gz > ../annotation/sgd_20191025_rna_coding.fasta

wget -nc -O ../annotation/Schizosaccharomyces_pombe_cds.fasta.gz https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/cds.fa.gz
gzip -dfc ../annotation/Schizosaccharomyces_pombe_cds.fasta.gz > ../annotation/Schizosaccharomyces_pombe_cds.fasta