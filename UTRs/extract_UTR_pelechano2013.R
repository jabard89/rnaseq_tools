# Created on 230421 by Jared Bard
# merges UTR annotations from yeasTSS and Pelechano_2013
library(tidyverse)
d_pelechano <- read_tsv("Pelechano2013_SuppData3_medianTranscripts.txt",comment="#") %>%
  mutate(UTR5_length = floor(median5),
         UTR3_length = floor(median3)) %>%
  rename('ORF'="gene") %>%
  select(ORF,UTR5_length,UTR3_length) %>% write_tsv("Pelechano2013_median_UTRs.tsv")
