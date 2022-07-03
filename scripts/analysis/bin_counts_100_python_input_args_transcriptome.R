library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
working.dir <- dirname(file)
name <- args[2]
print(paste0("Analyzing: ",file,"\n Exporting to: ",name,"_bin100.tsv.gz"))
# python numbers are 0-indexed, assuming the transcriptome came with 250nt UTR
d_raw <- read_tsv(file,comment="#") %>%
  mutate(Pos=Pos-250) %>%
  mutate(CDS_length = CDS_length-250-250)

d_normcounts <- d_raw %>%
  group_by(ORF) %>%
  mutate(counts.norm = Count/sum(Count,na.rm=T)) %>%
  ungroup %>%
  mutate(start_aligned = Pos) %>%
  select(ORF,start_aligned,counts.norm,CDS_length)

bin_orf <- function(df,N,vector) {
    bin_vector <- as.integer(cut_interval(vector,n=N,labels=seq(1,N)))
    bin_df <- tibble("start_aligned"=vector,"bin_vector"=bin_vector)
    df %>% full_join(bin_df,by="start_aligned")
}

bin_df <- function(df,bin) {
    if (bin=="UTR5") {
      df %>% mutate(bin_vector = -1) %>% return
    } else if (bin=="ORF") {
      df %>% group_by(ORF,CDS_length) %>%
        nest %>%
        mutate(data=map(data,bin_orf,100,seq(0,CDS_length-1))) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector-1) %>%
        ungroup() %>% return
    } else if (bin=="UTR3") {
      df %>% mutate(bin_vector = 100) %>% return
    } else {return(NULL)}
}

d_counts_binned <- d_normcounts %>%
    mutate(bin = case_when(start_aligned<0~"UTR5",
                           start_aligned>=0&start_aligned<CDS_length~"ORF",
                           start_aligned>=CDS_length~"UTR3")) %>%
    group_by(bin) %>%
    nest %>%
    mutate(data = case_when(bin=="UTR5"~map(data,bin_df,bin),
                            bin=="ORF"~map(data,bin_df,bin),
                            bin=="UTR3"~map(data,bin_df,bin),
                            TRUE ~ data)) %>%
    unnest(data) %>%
    ungroup %>%
    mutate(bin = bin_vector) %>%
    select(-bin_vector)
print("Writing")
d_binned_byORF <- d_counts_binned %>%
    group_by(ORF,bin) %>%
        summarise(counts = sum(counts.norm,na.rm=T)) %>%
    write_tsv(paste0(working.dir,"/",name,"_bin100.tsv.gz"))
