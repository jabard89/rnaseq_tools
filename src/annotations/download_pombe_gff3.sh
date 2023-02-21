#!/bin/bash

POMBE_GFF3="Schizosaccharomyces_pombe_all_chromosomes.gff3"
wget -nc -O "./${POMBE_GFF3}.gz" "https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"
gzip -dfc ./${POMBE_GFF3}.gz > ./${POMBE_GFF3}
