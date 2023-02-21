#!/bin/bash

POMBE_PREFIX_GFF3="Schizosaccharomyces_pombe_all_chromosomes_geneid_prefix.gff3"
SCER_GFF3="Saccharomyces_cerevisiae.R64-1-1.105_geneid.gff3"
OUTPUT_GFF3="Saccharomyces_cerevisiae.R64-1-1.105_merge_SPombe.gff3"

cat ./$SCER_GFF3 $POMBE_PREFIX_GFF3 > "./$(basename $SCER_GFF3 .gff3)_merge_$(basename $POMBE_PREFIX_GFF3)"