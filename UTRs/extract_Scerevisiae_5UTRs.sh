#!/bin/bash
curl -C - -o ./Saccharomyces_cerevisiae.R64-1-1.109.gff3.gz https://ftp.ensembl.org/pub/release-109/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gff3.gz
gzip -dfc ./Saccharomyces_cerevisiae.R64-1-1.109.gff3.gz > ./Saccharomyces_cerevisiae.R64-1-1.109.gff3
curl -C - -o ./Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gzip -dfc ./Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz > ./Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
curl -C - -o ./ScerYPDconsensusClusters.txt	http://yeastss.org/jbrowse/JBrowse_data/Saccharomyces_cerevisiae/ScerYPDconsensusClusters.txt
python extract-UTR5-yeasTSS.py Saccharomyces_cerevisiae.R64-1-1.109.gff3 Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeasTSS/ScerYPD.1.ctss yeasTSS/ScerYPD.2.ctss Saccharomyces_cerevisiae_YPD_YeasTSS_5UTRs.csv --organism=Scerevisiae --condition=YPD --max=900
Rscript --vanilla extract_UTR_pelechano2013.R
Rscript --vanilla merge_UTR_annotations.R
python gffutils_add_utrs.py --utr5_src=yeasTSS --utr3_src=Pelechano2013 Saccharomyces_cerevisiae.R64-1-1.109.gff3 Saccharomyces_cerevisiae_YPD_combined_yeasTSS_Pelechano2013_UTRs.tsv Saccharomyces_cerevisiae.R64-1-1.109_yeasTSS_Pelechano2013.gff3
python gffutils_add_utrs.py --utr5_src=Pelechano2013 --utr3_src=Pelechano2013 Saccharomyces_cerevisiae.R64-1-1.109.gff3 Pelechano2013_median_UTRs.tsv Saccharomyces_cerevisiae.R64-1-1.109_Pelechano2013.gff3
python ../scripts/annotation/extract_UTR_from_GFF.py Saccharomyces_cerevisiae.R64-1-1.109_Pelechano2013.gff3 Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa Saccharomyces_cerevisiae.R64-1-1.109_Pelechano2013_UTRs.csv --add_nts=20
