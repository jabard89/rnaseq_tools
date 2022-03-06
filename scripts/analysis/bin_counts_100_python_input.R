library(tidyverse)
# working.dir <- "/home/jbard/scratch-midway2/TSP"
# gene_labels <- read_tsv(paste0(working.dir,"/Analysis/src/201014_labeled_genes_scer.txt")) %>%
#   mutate(length.trans=length.prot*3)
# 
# files <- list.files(paste0(working.dir,"/big/hisat2/210614/counts"),
#                           pattern="*.tsv.gz",full.names = T)
gene_labels <- read_tsv("src/201014_labeled_genes_scer.txt") %>%
     mutate(length.trans=length.prot*3)
files <- list.files("../big/counts",pattern="*.tsv.gz",full.names=T)
names <- str_extract(files,pattern="(?<=_)S[0-9]*(?=_counts\\.tsv.gz)")
t <- tibble(names,files)

t$files %>% map2(t$names,function(file,name) {
  d_raw <- read_tsv(file) %>%
    left_join(gene_labels %>% select(ORF,length.trans),by="ORF") %>%
    rename("Counts"="Count") %>%
    filter(!is.na(length.trans))
  
  total.RPK <- d_raw %>% group_by(ORF,length.trans) %>%
      summarise(counts.perORF = sum(Counts,na.rm=T)) %>%
      ungroup() %>%
      mutate(RPK = counts.perORF / (length.trans/1000)) %>%
      summarise(totalRPK = sum(RPK,na.rm=T)) %>% pull(totalRPK)
  
  tpms <- d_raw %>% group_by(ORF,length.trans) %>%
      summarise(counts.perORF = sum(Counts,na.rm=T),.groups="drop") %>%
    ungroup %>%
      mutate(RPK = counts.perORF / (length.trans/1000),
             TPM = RPK*1e6/total.RPK)

  counts_per_ORF <- d_raw %>% group_by(ORF,length.trans) %>%
      summarise(counts.perORF = sum(Counts,na.rm=T),.groups="drop") %>%
      ungroup()
  
  d_normcounts <- d_raw %>%
    left_join(counts_per_ORF %>% select(ORF,counts.perORF),by="ORF") %>%
      mutate(counts.norm = Counts/counts.perORF) %>%
      mutate(start_aligned = Pos) %>%
      select(start_aligned,ORF,counts.norm)
  
  bin_orf <- function(df,N) {
    bin_vector <- as.integer(cut_interval(df$start_aligned,n=N,labels=seq(1,N)))
    df %>% add_column(bin_vector)
  }
  
  bin_df <- function(df,N,bin) {
    if (bin=="UTR5") {
      df %>% group_by(ORF) %>%
        nest %>%
        mutate(data=map(data,bin_orf,N)) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector-26) %>%
        ungroup() %>% return
    } else if (bin=="ORF") {
      df %>% group_by(ORF) %>%
        nest %>%
        mutate(data=map(data,bin_orf,N)) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector-1) %>%
        ungroup() %>% return
    } else if (bin=="UTR3") {
      df %>% group_by(ORF) %>%
        nest %>%
        mutate(data=map(data,bin_orf,N)) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector+99) %>%
        ungroup() %>% return
    } else {return(NULL)}
  }
  
  d_counts_binned <- d_normcounts %>%
    left_join(gene_labels %>% select(ORF,length.trans),by="ORF") %>%
    mutate(bin = case_when(start_aligned<0~"UTR5",
                           start_aligned>=0&start_aligned<length.trans~"ORF",
                           start_aligned>=length.trans~"UTR3")) %>%
    group_by(bin) %>%
    nest %>%
    mutate(data = case_when(bin=="UTR5"~map(data,bin_df,N=25,bin),
                            bin=="ORF"~map(data,bin_df,N=100,bin),
                            bin=="UTR3"~map(data,bin_df,N=25,bin),
                            TRUE ~ data)) %>%
    unnest(data) %>%
    ungroup %>%
    mutate(bin = bin_vector) %>%
    select(-bin_vector)

  d_binned_byORF <- d_counts_binned %>%
    group_by(ORF,bin) %>%
        summarise(counts = sum(counts.norm,na.rm=T)) %>%
    left_join(tpms %>% rename("TPM.hisat2" = "TPM") %>% select(ORF,TPM.hisat2),by=c("ORF")) %>%
    #write_tsv(paste0(working.dir,"/Analysis/output/bin100/counts_bin100_",name,".tsv.gz"))
    write_tsv(paste0("output/bin100/python/201706_TSP_counts_bin100_",name,".tsv.gz"))
})
