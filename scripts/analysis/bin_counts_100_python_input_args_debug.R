library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
args <- c("src/annotations/201014_labeled_genes_scer.txt","F:/Dropbox (Drummond Lab)/Data_JB/RNA-seq/JB/220301/star/HG49/HG49_pile.tsv",
          "HG49/","_pile.tsv")
file <- args[1]
working.dir <- dirname(file)
name <- str_extract(file,pattern=paste0("(?<=",args[2],").*(?=",args[3],")"))
gene_labels <- read_tsv("src/annotations/201014_labeled_genes_scer.txt") %>%
  mutate(length.CDS=length.prot*3)
print(paste0("Analyzing: ",name))
d_raw0 <- read_tsv(file,comment="#") %>%
    left_join(gene_labels %>% select(ORF,length.CDS),by="ORF") %>%
    rename("Counts"="Count") %>%
    filter(!is.na(length.CDS))

# find min and max
d_length <- d_raw0 %>% select(Pos,ORF,length.CDS) %>%
    group_by(ORF,length.CDS) %>% nest %>%
    mutate(start = map_dbl(data,function(df){min(df$Pos)}),
           end = map_dbl(data,function(df){max(df$Pos)})) %>%
    select(-data) %>% ungroup %>%
    mutate(start=if_else(start>-10,-10,start),
           end=if_else(end<(length.CDS+9),length.CDS+9,end),
           length.trans=end-start) #expand transcripts to at least the CDS length
d_raw <- d_raw0 %>%
    left_join(d_length %>% select(ORF,length.trans),by="ORF")
rm(d_raw0)
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

bin_orf <- function(df,N,vector) {
    bin_vector <- as.integer(cut_interval(vector,n=N,labels=seq(1,N)))
    bin_df <- tibble("start_aligned"=vector,"bin_vector"=bin_vector)
    df %>% left_join(bin_df,by="start_aligned")
}

bin_df <- function(df,bin) {
    if (bin=="UTR5") {
      df %>% group_by(ORF,start,end,length.CDS) %>%
        nest %>%
        mutate(data=map(data,bin_orf,10,seq(start,-1))) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector-11) %>%
        ungroup() %>% return
    } else if (bin=="ORF") {
      df %>% group_by(ORF,start,end,length.CDS) %>%
        nest %>%
        mutate(data=map(data,bin_orf,100,seq(0,length.CDS-1))) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector-1) %>%
        ungroup() %>% return
    } else if (bin=="UTR3") {
      df %>% group_by(ORF,start,end,length.CDS) %>%
        nest %>%
        mutate(data=map(data,bin_orf,10,seq(length.CDS,end))) %>%
        unnest(data) %>%
        mutate(bin_vector = bin_vector+99) %>%
        ungroup() %>% return
    } else {return(NULL)}
}
d_normcounts_test <- d_normcounts %>% filter(ORF=="YAL005C")
d_counts_binned <- d_normcounts_test %>%
    left_join(d_length,by="ORF") %>%
    mutate(bin = case_when(start_aligned<0~"UTR5",
                           start_aligned>=0&start_aligned<length.CDS~"ORF",
                           start_aligned>=length.CDS~"UTR3")) %>%
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
p1 <- ggplot(d_counts_binned%>%filter(ORF=="YAL005C"),
       aes(x=start_aligned,y=counts.norm))+
  geom_point()
p2 <- ggplot(d_binned_byORF%>%filter(ORF=="YAL005C"),
       aes(x=bin,y=counts))+
  geom_point()
library(cowplot)
plot_grid(p1,p2,ncol=1)
print("Writing")
d_binned_byORF <- d_counts_binned %>%
    group_by(ORF,bin) %>%
        summarise(counts = sum(counts.norm,na.rm=T)) #%>%
    left_join(tpms %>% rename("TPM.STAR" = "TPM") %>% select(ORF,TPM.STAR),by=c("ORF")) %>%
    left_join(tpms.trim %>% select(ORF,TPM.trim),by=c("ORF")) %>%
    write_tsv(paste0(working.dir,paste0("/",args[2],"bin100_",name,".tsv.gz")))
