setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/")

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
EUR_IPF <- fread("EUR_meta_IPF.b38.txt.gz")
Allen_ori <- fread("Allen_5way.txt.gz")
Allen <- Allen_ori %>% left_join(EUR_IPF, by=c("SNP"="MarkerName"))
Allen <- Allen %>% mutate(AF_meta = ifelse(Allele2.x == toupper(Allele1.y), Freq1, 1- Freq1 ))
Allen <- Allen %>% mutate(AF2 = case_when(AF_meta > 0.5 & AF_Allele2 > 0.5 ~ AF_Allele2,
                                       AF_meta <= 0.5 & AF_Allele2 <= 0.5 ~ AF_Allele2,
                                       AF_meta > 0.5 & AF_Allele2 <= 0.5 ~ 1 - AF_Allele2,
                                       AF_meta <= 0.5 & AF_Allele2 > 0.5 ~ 1 - AF_Allele2
))

Allen <- Allen %>% select(SNP, AF2)
Allen_ori <- Allen_ori %>% left_join(Allen, by="SNP")
Allen_ori <- Allen_ori %>% mutate(AF_Allele2 = AF2) %>% select(-AF2) %>% write.table("Allen_IPF.b38_EAFfiexed.txt", sep=" ", quote=F, row.names = F, col.names = T)

