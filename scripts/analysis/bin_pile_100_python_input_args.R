# Written 2022/04/06
# bins the counts within the CDS into 100 bins
# 221024 edited to only output CDS reads

library(tidyverse)
#args[1]: count_file (e.g. test1/test_pile.tsv.gz) should be .tsv.gz
#args[2]: sample_name
args <- commandArgs(trailingOnly = TRUE)
#args <- c("F:/Dropbox (Drummond Lab)/Data_JB/RNA-seq/JB/221019/snake/counts/JB125/221019_star_JB125_dist_test.tsv.gz",
#           "221019_star_JB125")
if (!file.exists(args[1])) {
    stop("file does not exists. Correct usage is bin_pile_100_python_input_args.R pile_file sample")
}
file <- args[1]
working.dir <- dirname(file)
sample <- args[2]
print(args)
if (is.na(sample)) {
  stop("sample not provided. Correct usage is bin_pile_100_python_input_args.R pile_file sample")
}
print(paste0("Analyzing: ",sample))

d_raw <- read_tsv(file,comment="#") %>%
  filter(Region=="CDS")

d_normcounts <- d_raw %>%
  group_by(ORF) %>%
  mutate(counts.norm = Count/sum(Count,na.rm=T)) %>%
  ungroup %>%
  mutate(start_aligned = Pos) %>%
  select(ORF,start_aligned,counts.norm,CDS_length)

bin_orf <- function(df,N,vector) {
  bin <- as.integer(cut_interval(vector,n=N,labels=seq(1,N)))
  bin_df <- tibble("start_aligned"=vector,"bin"=bin)
  df %>% full_join(bin_df,by="start_aligned")
}

d_counts_binned <- d_normcounts %>%
  group_by(ORF,CDS_length) %>%
  nest %>%
  mutate(data=map(data,bin_orf,100,seq(0,CDS_length-1))) %>%
  unnest(data) %>%
  ungroup
print("Writing")

d_binned_byORF <- d_counts_binned %>%
  group_by(ORF,bin) %>%
  summarise(counts = sum(counts.norm,na.rm=T)) %>%
  write_tsv(paste0(working.dir,"/",sample,"_bin100.tsv.gz"))
