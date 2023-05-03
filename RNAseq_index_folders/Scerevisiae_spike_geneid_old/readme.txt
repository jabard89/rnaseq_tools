# Jared Bard
# 230403
# creating annoations for mapping via STAR and kallisto
First ran bash download_Scerevisiae_genome_annotations.sh
That downloaded the latest genome release at the time (R64) from the Saccharomyces Genome Database (yeastgenome.org).
That combines the ORF coding transcripts, with the ncRNA sequences (tRNA, rRNA etc) and the spike-in transcript sequences.
I then ran gffutils_add_geneid.py to add geneids to every entry in the gff file (saccharomyces_cerevisiae_R64-3-1_20210421.geneid.gff).
I then copied the spike-in gff annotations into the newly annotated gff (saccharomyces_cerevisiae_R64-3-1_20210421.geneid.spike.gff).
Finally, I took the chromosome fasta files from the gff file and copied them into a separate file (saccharomyces_cerevisiae_R64-3-1_20210421_allChromosomes.fasta).
