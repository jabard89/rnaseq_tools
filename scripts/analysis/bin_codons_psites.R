#!/usr/bin/env Rscript --vanilla
# Jared Bard
# 4/18/2023
# Script to extract the reads from the beginning of ribosome profiling data

library(tidyverse)
library(optparse)

# Define the option parser
parser <- OptionParser()
parser = add_option(parser, c("-p", "--psites"), dest="psites", type="character", help="path to input psite file")
parser = add_option(parser, c("-m", "--minreads"), dest="read_threshold", type="integer", default=500, help="minimum read count threshold")
parser = add_option(parser, c("-f", "--firstn"), dest="first_n", type="integer", default=300, help="number of amino acids to output")
parser = add_option(parser, c("-b", "--binwidth"), dest="bin_width", type="integer", default=5, help="how many codons to bin")

# Parse the command line arguments
# args=c(paste0("--psites=/beagle3/dadrummond/jbard/riboseq_comp/riboWaltz/",
#               "psites/Iserman.2020_ribo_30C.control_W303_ribozero_BR2_psites.tsv.gz"),"--minreads=500","--firstn=300",
#        "--binwidth=5")
# options <- parse_args(parser,args)
options <- parse_args(parser)
print(options)
if(options$first_n %% options$bin_width != 0) {
  stop("Error: The bin width must divide evenly into the number of amino acids.")
}
basename <- str_extract(options$psites, "^.*(?=_psites.tsv.gz)")

# Read in the bedgraph file
col_specs <- list(
  "ORF" = col_character(),
  "end5" = col_double(),
  "psite" = col_double(),
  "end3" = col_double(),
  "length" = col_double(),
  "cds_start" = col_double(),
  "cds_stop" = col_double(),
  "psite_from_start" = col_double(),
  "psite_from_stop" = col_double(),
  "psite_region" = col_character()
)
psite_data <- vroom::vroom(options$psites,delim="\t",comment="#",
               col_names <- names(col_specs),skip=1,
               col_types = col_specs,
               col_select=c(ORF,psite,psite_region,cds_start,cds_stop)) %>%
    group_by(ORF,cds_start,cds_stop) %>% nest %>%
    mutate(data=map(data,function(df){
      df %>% filter(psite_region=="cds") %>%
        count(psite)
    })) %>%
    mutate(remove=map(data,function(df){
      if (sum(df$n)<options$read_threshold) {
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    })) %>%
    filter(remove==FALSE) %>% select(-remove)

# first group by codons, then extract the first n codons

d_codons <- psite_data %>%
  mutate(data=pmap(list(df=data,start=cds_start,stop=cds_stop),function(df,start,stop){
  	zero_tibble <- tibble(psite=seq(start,stop)) %>%
  	  mutate(codon=cut_width(psite,width=3,closed="right",center=start+1,labels=seq(1,(stop-start+1)/3)))%>%
  	  mutate(codon=as.integer(codon))
    df %>%
      full_join(zero_tibble,by=c("psite")) %>%
      replace_na(list(n=0)) %>%
      group_by(codon) %>%
      summarise(count=sum(n),.groups="drop") %>%
      mutate(count_norm=count/sum(count))
  })) %>%
  ungroup %>%
  select(-cds_start,-cds_stop) %>%
  unnest(c(data))
output_file1 <- paste0(basename,"_CDS_bycodon.tsv.gz")
vroom::vroom_write(d_codons,output_file1)

cut_tibble <- tibble(codon=seq(1,options$first_n)) %>%
  mutate(bin=cut_width(codon,width=options$bin_width,boundary=1,closed="left",
                       labels=seq(1,options$first_n/options$bin_width))) %>%
  mutate(bin=as.integer(bin))
d_front_downsample <- d_codons %>%
      filter(codon<=options$first_n) %>%
  left_join(cut_tibble,by="codon") %>%
  group_by(ORF,bin) %>%
  summarise(count=sum(count),
            count_norm=sum(count_norm),.groups='drop')
output_file2 <- paste0(basename, "_first",options$first_n,"_bin",options$bin_width,"_codon.tsv.gz")
vroom::vroom_write(d_front_downsample,output_file2)
