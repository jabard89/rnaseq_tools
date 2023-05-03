# Created on 230421 by Jared Bard
# merges UTR annotations from yeasTSS and Pelechano_2013
library(tidyverse)
d_pelechano <- read_tsv("Pelechano2013_SuppData3_medianTranscripts.txt",comment="#") %>%
  mutate(median5 = floor(median5),
         median3 = floor(median3)) %>%
  rename("ORF"="gene",
         "UTR3_length"="median3") %>%
  select(ORF,UTR3_length)
d_yeasTSS <- read_csv("Saccharomyces_cerevisiae_YPD_YeasTSS_5UTRs.csv",comment="#") %>%
  select(ORF,UTR5_length)

# d_combine <- full_join(d_pelechano,d_yeasTSS,by="ORF") %>%
#   mutate(UTR5_length=if_else(!is.na(UTR5_length_Pelechano),UTR5_length_Pelechano,UTR5_length_yeasTSS)) %>%
#   select(ORF,UTR5_length,UTR3_length) %>%
#   arrange(ORF)
d_combine <- full_join(d_pelechano,d_yeasTSS,by="ORF")
out_file <- "Saccharomyces_cerevisiae_YPD_combined_yeasTSS_Pelechano2013_UTRs.tsv"
header <- c("# Generated on 4/24/2023 by merge_UTR_annotations.R",
            "# combines 5' UTR annotations from yeasTSS (Saccharomyces_cerevisiae_YPD_YeasTSS_5UTRs.csv)",
            "# and 3' UTR annotations from Pelechano_2013 (Pelechano2013_SuppData3_medianTranscripts.txt)")
write_lines(header,out_file)
write_tsv(d_combine,out_file,append=T,col_names=T)
