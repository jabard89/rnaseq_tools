#!/bin/bash
fasta_name=Saccharomyces_cerevisiae.R64-1-1_transcriptome_padded_600up_30down_renamed
python ../../scripts/annotation/count_fasta_length.py ${fasta_name}.fasta ${fasta_name}.lengths.tsv
Rscript --vanilla add_utr_lengths.R ${fasta_name}.lengths.tsv 600 30 ${fasta_name}.utr_lengths.txt
