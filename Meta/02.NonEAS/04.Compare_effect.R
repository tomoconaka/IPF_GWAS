setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/")

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

# ALL <- fread("ALL_meta_IPF.b38.txt.gz")
# ALL <- ALL %>% mutate(missing_n = stringr::str_count(Direction, "\\?"))
# ALL <- ALL %>% filter(missing_n <= 5)
# ALL <- ALL %>% dplyr::select(-missing_n)
# 
# write.table(ALL, "ALL_meta_IPF.b38_filtered.txt", quote=F, col.names = T, row.names = F, sep="\t")
# ALL <- ALL %>% mutate(Pvalue = as.numeric(Pvalue))
# tmp <- ALL %>% filter(CHR == 1 & Pvalue < 5e-8)
# ALL %>% filter(CHR == 16 & POS >= 66900693 & Pvalue < 5e-8) %>% arrange(Pvalue) %>% head()
# ALL %>% filter(CHR == 16 & POS <= 1442406 + 500000 & POS > 1442406 - 500000  & Pvalue < 5e-8) %>% arrange(Pvalue) %>% head()
# ALL %>% filter(CHR == 2 & POS == 64652534)

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

ALL <- fread("ALL_meta_IPF.b38_filtered.txt.gz")
ALL <- ALL %>% filter(MarkerName %in% leadvariant)
ALL %>% write.table("significant_loci.txt", quote=F, col.names=T, sep="\t", row.names=F)
# EUR_IPF <- fread("EUR_meta_IPF.b38.txt.gz")
# EUR_IPF <- EUR_IPF %>% mutate(missing_n = stringr::str_count(Direction, "\\?"))
# EUR_IPF <- EUR_IPF %>% filter(missing_n <= 1)
# EUR_IPF <- EUR_IPF %>% dplyr::select(-missing_n)
# 
# write.table(EUR_IPF, "EUR_meta_IPF.b38_filtered.txt", quote=F, col.names = T, row.names = F, sep="\t")

EUR_IPF <- fread("EUR_meta_IPF.b38.txt.gz")

EUR_IPF <- EUR_IPF %>% mutate(Pvalue = as.numeric(Pvalue))
EUR_IPF <- EUR_IPF %>% filter(MarkerName %in% leadvariant)

EUR_IPF <- EUR_IPF %>% mutate(EA = ifelse(Effect > 0, toupper(Allele1), toupper(Allele2)),
                              NEA = ifelse(Effect > 0, toupper(Allele2), toupper(Allele1)),
                              EAF = ifelse(Effect > 0, Freq1, 1-Freq1),
                              Effect = ifelse(Effect > 0, Effect, Effect*(-1)),
                              Pvalue = Pvalue)

EUR_IPF <- EUR_IPF %>% dplyr::select(MarkerName, EA, NEA, EAF, Effect, StdErr, Pvalue)

# EAS_IPF <- fread("EAS_meta_IPF.b38.txt.gz")
# EAS_IPF <- EAS_IPF %>% mutate(missing_n = stringr::str_count(Direction, "\\?"))
# EAS_IPF <- EAS_IPF %>% filter(missing_n <= 0)
# EAS_IPF <- EAS_IPF %>% dplyr::select(-missing_n)
# 
# write.table(EAS_IPF, "EAS_meta_IPF.b38_filtered.txt", quote=F, col.names = T, row.names = F, sep="\t")

EAS_IPF <- fread("EAS_meta_IPF.b38.txt.gz")
EAS_IPF <- EAS_IPF %>% mutate(Pvalue = as.numeric(Pvalue))
EAS_IPF <- EAS_IPF %>% filter(MarkerName %in% leadvariant)

EAS_IPF <- EAS_IPF %>% left_join(EUR_IPF, by="MarkerName")
EAS_IPF <- EAS_IPF %>% mutate(Effect = ifelse(toupper(Allele1) == EA, Effect.x, Effect.x*(-1)),
                              StdErr = StdErr.x,
                              EAF = ifelse(toupper(Allele1) == EA, Freq1, 1 - Freq1),
                              Pvalue = Pvalue.x)
EAS_IPF <- EAS_IPF %>% mutate(Effect = ifelse(is.na(Effect), Effect.x, Effect),
                              EAF = ifelse(is.na(EAF), Freq1, EAF))

EAS_IPF <- EAS_IPF %>% dplyr::select(MarkerName, EAF, Effect, StdErr, Pvalue)


AFR_IPF <- fread("GBMI_AFR.txt.gz")
AFR_IPF <- AFR_IPF %>% filter(SNP %in% leadvariant)

AFR_IPF <- AFR_IPF %>% left_join(EUR_IPF, by=c("SNP"="MarkerName"))
AFR_IPF <- AFR_IPF %>% mutate(Effect = ifelse(toupper(Allele2) == EA, BETA, BETA*(-1)),
                              StdErr = SE,
                              EAF = ifelse(toupper(Allele2) == EA, AF_Allele2, 1 - AF_Allele2),
                              Pvalue = p.value,
                              MarkerName = SNP)

AFR_IPF <- AFR_IPF %>% dplyr::select(MarkerName, EAF, Effect, StdErr, Pvalue)


AMR_IPF <- fread("GBMI_AMR.txt.gz")
AMR_IPF <- AMR_IPF %>% filter(SNP %in% leadvariant)

AMR_IPF <- AMR_IPF %>% left_join(EUR_IPF, by=c("SNP"="MarkerName"))
AMR_IPF <- AMR_IPF %>% mutate(Effect = ifelse(toupper(Allele2) == EA, BETA, BETA*(-1)),
                              StdErr = SE,
                              EAF = ifelse(toupper(Allele2) == EA, AF_Allele2, 1 - AF_Allele2),
                              Pvalue = p.value,
                              MarkerName = SNP)

AMR_IPF <- AMR_IPF %>% dplyr::select(MarkerName, EAF, Effect, StdErr, Pvalue)

EUR_IPF <- EUR_IPF %>% dplyr::select(-EA, -NEA) %>% mutate(Ancestry = "EUR")
EAS_IPF <- EAS_IPF %>% mutate(Ancestry = "EAS")
AFR_IPF <- AFR_IPF %>% mutate(Ancestry = "AFR")
AMR_IPF <- AMR_IPF %>% mutate(Ancestry = "AMR")

ALL <- bind_rows(EUR_IPF, EAS_IPF, AFR_IPF, AMR_IPF)
ALL <- ALL %>% mutate(OR = exp(Effect),
                          LL = exp(Effect - qnorm(0.975)*StdErr),
                          UL = exp(Effect + qnorm(0.975)*StdErr))

ALL <- ALL %>% mutate(sig = ifelse(`Pvalue` < 0.05, FALSE, TRUE))
ALL <- ALL %>% mutate(CHR = as.numeric(str_split(MarkerName, pattern=":", simplify = T)[,1]),
                      POS = as.numeric(str_split(MarkerName, pattern=":", simplify = T)[,2]))
#data$CHR[data$rsid == "SLC39A11"] <- 17
#data$POS[data$rsid == "SLC39A11"] <- 72748046


labeli <- function(variable, value){
  names_li <- list("1:9109660" = "GPR157",
                   "1:150579566" = "MCL1*",
                   "1:155199564" = "MUC1/THSB3*",
                   "1:163372193" = "RGS5*",
                   "1:214484297" = "PTPN14",
                   "2:64652534" = "SERTAD2*",
                   "3:44804157" = "KIF15", 
                   "3:169774313" = "TERC", 
                   "4:88963935" = "FAM13A",
                   "5:1282299" = "TERT#", 
                   "5:169588475" = "SPDL1†", 
                   "6:7562999" = "DSP",
                   "6:27698141" = "HMGN4*",
                   "6:32436600" = "HLA-DRA*",
                   "6:43385242" = "ZNF318", 
                   "6:122431405" = "HSF2*",
                   "7:1936821" = "MAD1L1", 
                   "7:100032719" = "ZKSCAN1",
                   "8:119931204" = "DEPTOR", 
                   "9:106717987" = "ZNF462", 
                   "10:941334" = "GTPBP4*#",
                   "10:103920828" = "STN1", 
                   "11:1219991" = "MUC5B", 
                   "12:101892202" = "DRAM1*",
                   "13:82382479" = "chr13q31.1*§",
                   "13:112880670" = "ATP11A",
                   "15:40426335" = "IVD",
                   "15:85742593" = "AKAP13", 
                   "16:112241" = "NPRL3†", 
                   "16:67960834" = "SLC12A4§", 
                   "17:43370681" = "ARL4D*",
                   "17:46253848" = "KANSL1†", 
                   "17:76029575" = "EVPL*", 
                   "18:666625" = "TYMS*", 
                   "19:4717660" = "DPP9", 
                   "19:5847989" = "FUT3",
                   "20:63585923" = "RTEL1"
                   )
  return(names_li[value])
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
###Figure2B
ALL %>%
  mutate(name = fct_reorder(MarkerName, POS)) %>%
  mutate(name = fct_reorder(name, CHR)) %>%
  mutate(Ancestry = fct_relevel(Ancestry, "EAS", "EUR", "AFR", "AMR")) %>% 
  ggplot(aes(x=OR, xmin=LL, xmax=UL, y=Ancestry, color=Ancestry)) +
  facet_grid(~name, labeller = labeli) + scale_alpha_discrete(range = c(1, 0.3),guide = 'none') + 
  geom_point(aes(size=EAF, alpha=sig)) + theme_minimal()  +
  geom_errorbarh(aes(height = 0, alpha=sig)) +  xlab("Odds ratio") + scale_x_log10() + 
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15), 
        strip.text.x = element_text(size = 13, angle = 90, face="italic"),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom") +
  coord_flip() + geom_vline(xintercept=1, alpha=0.5) +
  labs(size="Effect Allele Frequency", color = "Ancestry") + 
  scale_color_manual(labels = c("EAS", "EUR", "AFR", "AMR"), values=c("#e4007f", "#0068b7", "#f39800", "#009944"))
ggsave(paste0("../repo/IPF_GWAS/Figures/Fig1c.png"), width=50, height=10, units = "cm", dpi=300)




