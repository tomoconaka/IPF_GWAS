setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/")

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

leadvariant <- c("1:9109660", "1:150579566", "1:155199564", "1:163372193","1:214484297", 
                 "2:64652534",
                 "3:44804157", "3:169774313", 
                 "4:88963935",
                 "5:1282299", "5:169588475", 
                 "6:7562999", "6:27698141","6:32436600","6:43385242","6:122431405", 
                 "7:1936821", "7:100032719",
                 "8:119931204", 
                 "9:106717987", 
                 "10:941334","10:103920828", 
                 "11:1219991",
                 "12:101892202",
                 "13:82382479", "13:112880670",
                 "15:40426335", "15:85742593", 
                 "16:112241", "16:67960834", 
                 "17:43370681", "17:46253848", "17:76029575", 
                 "18:666625", 
                 "19:4717660", "19:5847989",
                 "20:63585923"
)

ALL_IPF <- fread("ALL_meta_IPF.b38_filtered.txt.gz")
ALL_IPF <- ALL_IPF %>% mutate(Pvalue = as.numeric(Pvalue))
ALL_IPF <- ALL_IPF %>% dplyr::filter(MarkerName %in% leadvariant)

EUR_IPF <- fread("EUR_meta_IPF.b38.txt.gz")
EUR_IPF <- EUR_IPF %>% mutate(Pvalue = as.numeric(Pvalue))
EUR_IPF <- EUR_IPF %>% dplyr::filter(MarkerName %in% leadvariant)

EAS_IPF <- fread("EAS_meta_IPF.b38.txt.gz")
EAS_IPF <- EAS_IPF %>% mutate(Pvalue = as.numeric(Pvalue))
EAS_IPF <- EAS_IPF %>% dplyr::filter(MarkerName %in% leadvariant)

AFR_IPF <- fread("GBMI_AFR.txt.gz")
AFR_IPF <- AFR_IPF %>% dplyr::filter(SNP %in% leadvariant)

AMR_IPF <- fread("GBMI_AMR.txt.gz")
AMR_IPF <- AMR_IPF %>% dplyr::filter(SNP %in% leadvariant)

ALL_IPF <- ALL_IPF %>% mutate(EA = ifelse(Effect > 0, toupper(Allele1), toupper(Allele2)),
                              NEA = ifelse(Effect > 0, toupper(Allele2), toupper(Allele1)),
                              Effect = ifelse(Effect > 0, Effect, Effect*(-1)),
                              Pvalue.ALL = Pvalue,
                              N.ALL = TotalSampleSize)

ALL_IPF <- ALL_IPF %>% mutate(OR.ALL = paste0(round(exp(Effect), 2), " [", round(exp(Effect + qnorm(0.025)*StdErr), 2),": ",round(exp(Effect + qnorm(0.975)*StdErr), 2),"]"))

ALL_IPF <- ALL_IPF %>% dplyr::select(CHR, POS, MarkerName, EA, NEA, OR.ALL, Pvalue.ALL, N.ALL)

EAS_IPF <- EAS_IPF %>% left_join(ALL_IPF, by=c("MarkerName"))

EAS_IPF <- EAS_IPF %>% mutate(Effect = ifelse(toupper(Allele1) == EA, Effect, Effect*(-1)),
                              Pvalue.EAS = Pvalue,
                              N.EAS = TotalSampleSize,
                              EAF.EAS = ifelse(toupper(Allele1) == EA, Freq1, 1- Freq1)
) %>% 
  mutate(OR.EAS = paste0(round(exp(Effect), 2), " [", round(exp(Effect + qnorm(0.025)*StdErr), 2),": ",round(exp(Effect + qnorm(0.975)*StdErr), 2),"]")) %>% 
  dplyr::select(MarkerName, OR.EAS, Pvalue.EAS, EAF.EAS, N.EAS)

ALL <- left_join(ALL_IPF, EAS_IPF, by="MarkerName")

EUR_IPF <- EUR_IPF %>% left_join(ALL_IPF, by=c("MarkerName"))

EUR_IPF <- EUR_IPF %>% mutate(Effect = ifelse(toupper(Allele1) == EA, Effect, Effect*(-1)),
                              Pvalue.EUR = Pvalue,
                              N.EUR = TotalSampleSize,
                              EAF.EUR = ifelse(toupper(Allele1) == EA, Freq1, 1- Freq1)
                              ) %>% 
  mutate(OR.EUR = paste0(round(exp(Effect), 2), " [", round(exp(Effect + qnorm(0.025)*StdErr), 2),": ",round(exp(Effect + qnorm(0.975)*StdErr), 2),"]")) %>% 
  dplyr::select(MarkerName, OR.EUR, Pvalue.EUR, EAF.EUR, N.EUR)

ALL <- left_join(ALL, EUR_IPF, by="MarkerName")

AFR_IPF <- AFR_IPF %>% left_join(ALL_IPF, by=c("SNP"="MarkerName"))

AFR_IPF <- AFR_IPF %>% mutate(Effect = ifelse(toupper(Allele2) == EA, BETA, BETA*(-1)),
                              Pvalue.AFR = p.value,
                              N.AFR = N,
                              EAF.AFR = ifelse(toupper(Allele2) == EA, AF_Allele2, 1-AF_Allele2)
) %>% 
  mutate(OR.AFR = paste0(round(exp(Effect), 2), " [", round(exp(Effect + qnorm(0.025)*SE), 2),": ",round(exp(Effect + qnorm(0.975)*SE), 2),"]")) %>% 
  dplyr::select(SNP, OR.AFR, Pvalue.AFR, EAF.AFR, N.AFR)

ALL <- left_join(ALL, AFR_IPF, by=c("MarkerName"="SNP"))


AMR_IPF <- AMR_IPF %>% left_join(ALL_IPF, by=c("SNP"="MarkerName"))

AMR_IPF <- AMR_IPF %>% mutate(Effect = ifelse(toupper(Allele2) == EA, BETA, BETA*(-1)),
                              Pvalue.AMR = p.value,
                              N.AMR = N,
                              EAF.AMR = ifelse(toupper(Allele2) == EA, AF_Allele2, 1-AF_Allele2)
) %>% 
  mutate(OR.AMR = paste0(round(exp(Effect), 2), " [", round(exp(Effect + qnorm(0.025)*SE), 2),": ",round(exp(Effect + qnorm(0.975)*SE), 2),"]")) %>% 
  dplyr::select(SNP, OR.AMR, Pvalue.AMR, EAF.AMR, N.AMR)

ALL <- left_join(ALL, AMR_IPF, by=c("MarkerName"="SNP"))


ALL_IPF <- fread("ALL_meta_IPF.percontinental.b38.txt.gz")
ALL_IPF <- ALL_IPF %>% mutate(Pvalue = as.numeric(Pvalue))
ALL_IPF <- ALL_IPF %>% dplyr::filter(MarkerName %in% leadvariant)

ALL_IPF <- ALL_IPF %>% dplyr::select(MarkerName, HetPVal)

ALL <- ALL %>% left_join(ALL_IPF, by="MarkerName")

write.table(ALL, "../repo/IPF_GWAS/Tables/Table3_ALL_meta_IPF.tsv", sep="\t", quote=F, row.names = F, col.names = T)



