#!/bin/bash
SAMPLE_INDEX=$1
PROJECTFOLDER="/home/jbard/beagle3-dadrummond/jbard/211216"
star_bin="/home/jbard/tools/STAR-2.7.10a/bin/Linux_x86_64_static"
samtools_bin="/home/jbard/tools/samtools-1.15"
OUTPUT="$PROJECTFOLDER/STAR/output_220303"
INDEX_IN="$PROJECTFOLDER/STAR/index_genome/input"
INDEX_OUT="$PROJECTFOLDER/STAR/index_genome/merged"

$star_bin/STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir $INDEX_OUT \
	--genomeFastaFiles $INDEX_IN/Scer.R64_Spombe.ASM294v2_merge.fa \
	--sjdbGTFfile $INDEX_IN/Scer.R64_Spombe.ASM294v2_merge.gff3 \
	--sjdbOverhang 99 \
	--genomeSAindexNbases 10
