# Psites to bedgraph
# script to convert the output of ribowaltz p-sites to bedgraph format
# Jared Bard
# 230404
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE) # psite file, length of UTR padding, and output files
psite_file <- args[1]
utr_length <- as.integer(args[2])
out_all_file <- args[3]
out_inframe_file <- args[4]

df_input <- read_tsv(psite_file) %>%
  mutate(psite=psite-1) %>% # make it 0-based
  mutate(in_frame = if_else((psite-utr_length)%%3==0,TRUE,FALSE))

writeLines(c(paste0("# P-site footprints for all reads"),
             paste0("# Generated ",Sys.time()," by riboWaltz_psite_to_bedgraph"),
             paste0("# Input = ",psite_file,", UTR_padding=",utr_length),
             "track type=bedGraph name=footprints graphType=bar visibility=full autoScale=on"),
           sep="\n",
           con=out_all_file)
out_all <- df_input %>% 
  group_by(transcript) %>%
  count(psite) %>%
  mutate(psite2=psite) %>%
  select(transcript,psite,psite2,n) %>%
  write_tsv(out_all_file,append=T,col_names=F)

writeLines(c(paste0("# P-site footprints just for in-frame reads"),
             paste0("# Generated ",Sys.time()," by riboWaltz_psite_to_bedgraph"),
             paste0("# Input = ",psite_file,", UTR_padding=",utr_length),
             "track type=bedGraph name=footprints graphType=bar visibility=full autoScale=on"),
           sep="\n",
           con=out_inframe_file)
out_inframe <- df_input %>%
  filter(in_frame==TRUE) %>%
  group_by(transcript) %>%
  count(psite) %>%
  mutate(psite2=psite) %>%
  select(transcript,psite,psite2,n) %>%
  write_tsv(out_inframe_file,append=T,col_names=F)


