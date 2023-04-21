library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)
df <- read_tsv(args[1],comment="#") %>%
mutate(UTR5_length=as.integer(args[2]),UTR3_length=as.integer(args[3])) %>%
mutate(CDS_length=Length-UTR5_length-UTR3_length) %>%
select(ORF,UTR5_length,CDS_length,UTR3_length) %>%
arrange(ORF) %>%
write_tsv(args[4],col_names=FALSE)
