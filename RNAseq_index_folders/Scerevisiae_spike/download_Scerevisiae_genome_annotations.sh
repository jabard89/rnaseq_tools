#!/bin/bash
DIR=.
curl -C - -o $DIR/S288C_reference_genome_R64-3-1_20210421.tgz http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-3-1_20210421.tgz
tar -xvzf $DIR/S288C_reference_genome_R64-3-1_20210421.tgz

gzip -dfc $DIR/S288C_reference_genome_R64-3-1_20210421/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz > $DIR/saccharomyces_cerevisiae_R64-3-1_20210421.gff
gzip -dfc $DIR/S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421.fasta.gz > $DIR/orf_coding_all_R64-3-1_20210421.fasta
gzip -dfc $DIR/S288C_reference_genome_R64-3-1_20210421/rna_coding_R64-3-1_20210421.fasta.gz > $DIR/rna_coding_R64-3-1_20210421.fasta

cat spike.fasta $DIR/orf_coding_all_R64-3-1_20210421.fasta $DIR/rna_coding_R64-3-1_20210421.fasta > saccharomyces_cerevisiae_R64-3-1_2021041_orf_rna_spike.fasta






