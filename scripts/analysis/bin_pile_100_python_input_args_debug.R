library(tidyverse)
#args[1]: gene_labels (e.g. src/201014_gene_labels.tsv)
#args[2]: count_file (e.g. test1/test_pile.tsv.gz)
#args[3]: sample_prefix (if no prefix, use "")
#args[4]: sample_suffix (e.g. .tsv.gz)
args <- commandArgs(trailingOnly = TRUE)
# args <- c("src/annotations/201014_labeled_genes_scer.txt",
#           "F:/Dropbox (Drummond Lab)/Data_JB/RNA-seq/JB/220301/star/HG53/220403_STAR_HG53_pile_test.tsv.gz",
#           "220403_STAR_","_pile_test.tsv.gz")
for (f in args[1:2]) {
  if (!file.exists(f)) {
    stop("file does not exists. Correct usage is bin_pile_100_python_input_args.R gene_labels pile_file sample_prefix sample_suffix")
  }
}
file <- args[2]
working.dir <- dirname(file)
if (args[3]=="") {
  name <- str_extract(file,pattern=paste0("[^/]*(?=",args[4],")"))
} else {
  name <- str_extract(file,pattern=paste0("(?<=",args[3],").*(?=",args[4],")"))
}
print(args)
if (is.na(name)) {
  stop("Name not extracted properly")
}
gene_labels <- read_tsv(args[1]) %>%
  mutate(length.CDS=length.prot*3)
print(paste0("Analyzing: ",name))
d_raw0 <- read_tsv(file,comment="#") %>%
  filter(Pos>=0) %>%
    left_join(gene_labels %>% select(ORF,length.CDS),by="ORF") %>%
    rename("Counts"="Count")
# find min and max
d_length <- d_raw0 %>% select(Pos,ORF,length.CDS) %>%
    group_by(ORF,length.CDS) %>% nest %>% ungroup %>%
    mutate(end = map_dbl(data,function(df){max(df$Pos)})) %>%
    mutate(length.CDS = if_else(is.na(length.CDS),end,length.CDS)) %>%
    select(ORF,length.CDS) %>% ungroup

d_raw <- d_raw0 %>%
  select(-length.CDS) %>%
  left_join(d_length %>% select(ORF,length.CDS),by="ORF") %>%
  filter(Pos<length.CDS)
rm(d_raw0)
counts_per_ORF <- d_raw %>% group_by(ORF) %>%
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
  print(bin)
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

d_counts_binned <- d_normcounts %>%
    left_join(d_length,by="ORF") %>%
    group_by(ORF,length.CDS) %>%
    nest %>%
    mutate(data = map(data,bin_orf,100,seq(0,length.CDS-1)))%>%
    unnest(data) %>%
    ungroup %>%
    mutate(bin = bin_vector) %>%
    select(-bin_vector)

print("Writing")
d_binned_byORF <- d_counts_binned %>%
    group_by(ORF,bin) %>%
        summarise(counts = sum(counts.norm,na.rm=T))%>%
    write_tsv(paste0(working.dir,paste0("/",args[3],"bin100_",name,".tsv.gz")))
