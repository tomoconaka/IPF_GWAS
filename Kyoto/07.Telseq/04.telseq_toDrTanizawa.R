setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/")
library(data.table)
library(tidyverse)

telseq <- fread("Telseq/Telseq.summary.gz", sep=" ")
ipf <- fread("../data/IPsample20220601.list")
ipf <- ipf %>% select(`CGM ID (二次)`, WGS)
colnames(ipf) <- c("ID", "WGS_platform")

telseq <- telseq %>% left_join(ipf, by="ID")
telseq <- telseq %>% mutate(WGS_platform = case_when(ID %in% c("PFKT2312_wd", "PFKT2334_wd", "PFKT2377_wd") ~ "DNBSeq G400RS",
                                                     is.na(WGS_platform) ~  "HiSeq X 15x",
                                                     TRUE ~ WGS_platform))

telseq <- telseq %>% filter(grepl("^TL|^PF", ID)) %>% filter(WGS_platform %in% c("HiSeq X 15x", "NovaSeq 6000 15x"))

telseq %>% openxlsx::write.xlsx("Telseq/Telseq_toDrTanizawa.xlsx")

