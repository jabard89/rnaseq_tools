#!/bin/bash
SAMPLE_INDEX=$1
PROJECTFOLDER="/home/jbard/beagle3-dadrummond/jbard/211216"
star_bin="/home/jbard/tools/STAR-2.7.10a/bin/Linux_x86_64_static"
OUTPUT="$PROJECTFOLDER/STAR/output_220303"
INDEX="$PROJECTFOLDER/STAR/index_genome/merged"
SAMPLE_FASTQS="$PROJECTFOLDER/fastq/220103_A00639_0973_AHTCJMDRXY-ADr-JB-1S-RS/FastQ"
NTHREADS=8
SNAME="211216_STAR_JB${SAMPLE_INDEX}"
INPUT_R1="${SAMPLE_FASTQS}/ADr-JB0${SAMPLE_INDEX}*_R1_001.fastq.gz"
INPUT_R2="${SAMPLE_FASTQS}/ADr-JB0${SAMPLE_INDEX}*_R2_001.fastq.gz"

echo "SAMPLE_INDEX=${SAMPLE_INDEX}";ls $INPUT_R1;ls $INPUT_R2

$star_bin/STAR \
	--runThreadN $NTHREADS \
	--genomeDir $INDEX \
	--readFilesCommand zcat \
	--readFilesIn $INPUT_R1 $INPUT_R2 \
	--outFileNamePrefix "${OUTPUT}/JB${SAMPLE_INDEX}/${SNAME}_" \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN $NTHREADS \
	--alignMatesGapMax 20000 \
	--limitBAMsortRAM 1445804817
