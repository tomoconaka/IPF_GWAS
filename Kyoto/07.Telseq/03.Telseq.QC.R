setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/")

library(data.table)
library(tidyverse)

ipf <- fread("../data/IPsample20220601.list")
telseq <- fread("Telseq/Telseq.summary.gz", sep=" ")
ipf <- ipf %>% select(`CGM ID (二次)`, WGS)
colnames(ipf) <- c("ID", "WGS_platform")

telseq <- telseq %>% left_join(ipf, by="ID")
telseq <- telseq %>% mutate(WGS_platform = case_when(ID %in% c("PFKT2312_wd", "PFKT2334_wd", "PFKT2377_wd") ~ "DNBSeq G400RS",
                                                     is.na(WGS_platform) ~  "HiSeq X 15x",
                                                     TRUE ~ WGS_platform))

pheno <- fread("pheno.txt")
telseq <- telseq %>% inner_join(pheno, by=c("ID"="FID"))
telseq <- telseq %>% mutate(DNB = ifelse(WGS_platform %in% c("DNBSeq G400RS", "DNBSeq-T7 15x"), TRUE, FALSE))

png("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/SupFig4.png", width = 600)
ggplot(telseq, aes(x=DNB, y=telomere)) +
  geom_jitter(aes(color=DNB)) + geom_violin(aes(fill=DNB), trim = FALSE) + geom_boxplot(alpha=0.6,trim = FALSE) +
  scale_fill_manual(values=c("#00AFBB", "#E7B800"), name = "Sequencing platforms", labels = c(paste0("Illumina (N=",table(telseq$DNB)[1],")"), paste0("MGI (N=",table(telseq$DNB)[2],")"))) +
  scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme_classic() + ylab("Estimated telomere length (kb)") +
  theme(legend.title = element_text( size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(size=20)) +
  guides(color = FALSE)
dev.off()

ipf <- fread("../data/IPsample20220601.list")
ipf <- ipf %>% select(`CGM ID (二次)`, age)
colnames(ipf) <- c("ID", "age")
ipf <- ipf %>% mutate(age = as.numeric(age))
agp3000 <- fread("AGP3K_agesexsmoking20210819.txt")
agp3000 <- agp3000 %>% select(ID, age)
agp3000 <- agp3000 %>% mutate(age = as.numeric(age))
tmp <- bind_rows(ipf, agp3000)
tmp1 <- telseq %>% left_join(tmp, by="ID")
tmp1 <- tmp1 %>% drop_na(age.y)
tmp1 <- tmp1 %>% filter(DNB == FALSE)

LM <- glm(telomere ~ age.y*geneticSex, data=tmp1, family = "gaussian")

png("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/SupFig3.png", width = 700)
ggplot(tmp1, aes(x=age.y, y=telomere, col=geneticSex)) + geom_jitter() + xlab("age") +
  theme_classic() + ylab("Estimated telomere length (kb)") + 
  geom_smooth(method=lm, se=TRUE, mapping = aes(colour = geneticSex)) +
  ggtitle(paste0("N=",dim(tmp1)[1])) +
  theme(legend.title = element_text( size = 20), 
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  annotate("text", x = 75, y = 10.5, size=5, 
           label = paste0("Female: ",round(summary(LM)$coefficients[2,1], 3),"kb/yr\nMale: ",round(summary(LM)$coefficients[2,1] + summary(LM)$coefficients[4,1], 3),"kb/yr"))
dev.off()

LM <- glm(telomere ~ age.y*geneticSex, data=tmp1, family = "gaussian")
tmp1$tel_sd <- scale(residuals(LM))

ggplot(tmp1, aes(x=as.factor(ILD), y=tel_sd, color=as.factor(ILD))) + geom_violin(trim = FALSE) +
  geom_jitter() + theme_classic() + ylab("telomere length") +
  geom_boxplot() + theme(legend.title = element_text( size = 20),
                         legend.text = element_text(size = 20),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20)) 

tmp1 %>% write.table("telomere.txt", quote=F, col.names = T, row.names = F, sep="\t")
