#!/usr/bin/env Rscript
# Jared Bard
# 230312
# based on code found here: https://www.biostars.org/p/83901/
# gtf needs to have an exon feature for every gene, which can be added using the script gffutils_fix_missing_exon.py
# usage: pull_gene_lengths_from_gtf.R -i input.gtf -o output.tsv

# load libraries

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GenomicFeatures"))

option_list <- list( 
    make_option(c("-i", "--input_gtf"), type="character", default=NULL,help="input gtf file"),
	make_option(c("-o", "--output"), type="character",default=NULL,help="output tsv file")
)
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input_gtf)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file).n", call.=FALSE)
}

# write out header
fileConn<-file(opt$output)
writeLines(c("# Output of pull_gene_lengths_from_gtf.R",
			paste0("# input gtf was ",opt$input_gtf),
			paste0("# file written ",Sys.time())), fileConn)
close(fileConn)
# read in gtf file

txdb <- makeTxDbFromGFF(opt$input_gtf,format="gtf")
CDS.list.per.gene <- cdsBy(txdb,by="gene")
CDS.gene.sizes <- sum(width(GenomicRanges::reduce(CDS.list.per.gene)))
gene_lengths <- tibble("gene_id"=names(CDS.gene.sizes),"merged.exon.length"=CDS.gene.sizes)
write_tsv(gene_lengths,opt$output,append=TRUE,col_names=TRUE)