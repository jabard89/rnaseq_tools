# Run riboWaltz
# Jared Bard
# 07/12/2022
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE) # input sample name
sample <- args[1]
d_sample <- tibble("Sample"=sample) %>%
  separate(Sample,into=c("Paper","Type","Condition","Strain","Depletion","Biorep"),
           sep="_",remove=F)
stopifnot(d_sample$Type[1]=="ribo")

src.dir <- getwd()

d_lengths <- read_tsv("~/repos/rnaseq_tools/src/annotations/200325-scer-features.txt",
                      comment="#") %>% filter(classification=="Verified")
annot_dt <- data.table::data.table(transcript=d_lengths$ORF,l_tr=d_lengths$length.nt+500,
                                   l_utr5=250,l_cds=d_lengths$length.nt,l_utr3=250)

bam_file <- paste0(sample,"_transcriptome_Aligned.sortedByCoord.out")
sample_folder <- paste0(src.dir,"/mapped_reads/",sample)
names(sample) <- bam_file

temp_list <- riboWaltz::bamtolist(bamfolder = sample_folder,
                                  annotation=annot_dt,
                                  transcript_align=TRUE,
                                  name_samples=sample)

heatmap <- riboWaltz::rends_heat(temp_list,annot_dt,
                                 "sample"=sample)
pdf(paste0(src.dir,"/riboWaltz/ribogrid/",sample,"_ribogrid.pdf"))
heatmap[["plot"]]+labs(title=sample)+
  theme(plot.title=element_text(size=12))+
  scale_fill_viridis_c()
dev.off()

offset <- riboWaltz::psite(temp_list)
offset_file <- paste0(src.dir,"/riboWaltz/offset/",sample,"_offset.tsv")
writeLines(c(paste0("# P site offsets"),
             paste0("# Generated ",Sys.time()," by Jared Bard"),
             paste0("# Analyzed by riboWaltz_count.R using riboWaltz=",
                    packageVersion("riboWaltz")),
             paste0("# Input = ",sample_folder,"/",bam_file)),
           sep="\n",
           con=offset_file)
offset %>% as_tibble %>% write_tsv(offset_file,append=TRUE,col_names=TRUE)

psite_updated <- riboWaltz::psite_info(temp_list,offset)
rm(temp_list)
data.table::fwrite(psite_updated[[1]],
                   file=paste0(src.dir,"/riboWaltz/psites/",sample,"_psites.tsv.gz"),
                   sep="\t",row.names=FALSE,col.names=TRUE,compress="auto")

full_coverage <- riboWaltz::cds_coverage(psite_updated,annot_dt)
full_coverage_file <- paste0(src.dir,"/riboWaltz/coverage/",sample,"_full_count.tsv")
writeLines(c(paste0("# Ribosome coverage"),
             paste0("# Generated ",Sys.time()," by Jared Bard"),
             paste0("# Analyzed by riboWaltz_count.R using riboWaltz=",
                    packageVersion("riboWaltz")," cds_coverage with no trimming"),
             paste0("# Input = ",sample_folder,"/",bam_file)),
           sep="\n",
           con=full_coverage_file)
d_full_coverage <- full_coverage %>% as_tibble %>%
  rename("ORF"="transcript")
names(d_full_coverage)[names(d_full_coverage)==sample] <- "Count"
d_full_coverage <- d_full_coverage %>%
  mutate(FPKB = Count/(length_cds/1000)) %>%
  mutate(TPM = FPKB/(sum(FPKB)/1e6)) %>%
  select(ORF,Count,TPM) %>%
  write_tsv(full_coverage_file,append=TRUE,col_names=TRUE)

trimmed_coverage <- riboWaltz::cds_coverage(psite_updated,annot_dt,
                                            start_nts=100,stop_nts=30)
trimmed_coverage_file <- paste0(src.dir,"/riboWaltz/coverage/",sample,
                                "_trim_100nt_to_-30nt_count.tsv")
writeLines(c(paste0("# Ribosome coverage"),
             paste0("# Generated ",Sys.time()," by Jared Bard"),
             paste0("# Analyzed by riboWaltz_count.R using riboWaltz=",
                    packageVersion("riboWaltz")," cds_coverage with start_nts=100, stop_nts=30"),
             paste0("# Input = ",sample_folder,"/",bam_file)),
           sep="\n",
           con=trimmed_coverage_file)
d_trimmed_coverage <- trimmed_coverage %>% as_tibble %>%
  rename("ORF"="transcript")
names(d_trimmed_coverage)[names(d_trimmed_coverage)==sample] <- "Count"
d_trimmed_coverage <- d_trimmed_coverage %>%
  mutate(FPKB = Count/(length_selection/1000)) %>%
  mutate(TPM = FPKB/(sum(FPKB)/1e6)) %>%
  select(ORF,Count,TPM) %>%
  write_tsv(trimmed_coverage_file,append=TRUE,col_names=TRUE)

metaprofile <- riboWaltz::metaprofile_psite(psite_updated,annot_dt,
                                            "sample"=sample)
pdf(paste0(src.dir,"/riboWaltz/metaprofile/",sample,"_metaprofile.pdf"))
metaprofile[[paste0("plot_",sample)]]+labs(title=sample)+
  theme(plot.title=element_text(size=12))
dev.off()


