setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/")

library(data.table)
library(tidyverse)
agp3000 <- fread("../data/AGP3000sample.txt")

tmp1 <- fread("fail.verifyBamID2.sample")

agp3000 <- agp3000 %>% mutate(VerifyBamID = ifelse(FID %in% tmp1$V1, TRUE, FALSE))

tmp2 <- fread("fail.mappingrate0.97.chimera0.08.sample")
agp3000 <- agp3000 %>% mutate(MappingRate0.97 = ifelse(FID %in% tmp2$V1[tmp2$V2 < 0.97], TRUE, FALSE),
                              Chimera0.08 = ifelse(FID %in% tmp2$V1[tmp2$V3 > 0.08], TRUE, FALSE))

tmp3 <- fread("fail.XXY.sample")
agp3000 <- agp3000 %>% mutate(XXY = ifelse(FID %in% tmp3$V1, TRUE, FALSE))

agp3000 <- agp3000 %>% mutate(group = "Control")
agp3000 <- agp3000 %>% mutate(sex = ifelse(sex == 2, "Female", "Male"))

agp3000 <- agp3000 %>% mutate(tokeep = ifelse(VerifyBamID == FALSE & MappingRate0.97 == FALSE & Chimera0.08 == FALSE & XXY == FALSE, 1, 0))
agp3000 <- agp3000 %>% mutate(WGS = "HiSeq X 15x")
#agp3000 <- agp3000 %>% filter(tokeep == 1)

ipf <- fread("../data/IPsample20220601.list")
ipf <- ipf %>% mutate(FID = `CGM ID (二次)`,
                      IID = `CGM ID (二次)`) %>% select(-`CGM ID (二次)`)

ipf <- ipf %>% mutate(VerifyBamID = ifelse(FID %in% tmp1$V1, TRUE, FALSE))
ipf <- ipf %>% mutate(MappingRate0.97 = ifelse(FID %in% tmp2$V1[tmp2$V2 < 0.97], TRUE, FALSE),
                              Chimera0.08 = ifelse(FID %in% tmp2$V1[tmp2$V3 > 0.08], TRUE, FALSE))
ipf <- ipf %>% mutate(XXY = ifelse(FID %in% tmp3$V1, TRUE, FALSE))
ipf <- ipf %>% mutate(tokeep = ifelse(VerifyBamID == FALSE & MappingRate0.97 == FALSE & Chimera0.08 == FALSE & XXY == FALSE, 1, 0))


ipf <- ipf %>% mutate(Project = ifelse(grepl("^TL", FID), "TL", "PF"),
                      Facility = faciliy,
                      sex =  ifelse(sex == 1, "Female", "Male"),
                      group = case_when(grepl("1", `diagnostic code`) & grepl("2", `diagnostic code`) ~ "fIPF",
                                        grepl("2", `diagnostic code`) ~ "FPF",
                                        grepl("0", `diagnostic code`) ~ "Control",
                                        grepl("1", `diagnostic code`) ~ "sIPF")
                      ) 

#ipf <- ipf %>% filter(tokeep == 1)
table(ipf[,c("WGS", "group")])

ipf1 <- ipf %>% select(any_of(colnames(agp3000)))

all <- bind_rows(agp3000, ipf1)
#all <- all %>% filter(tokeep == 1)
all <- all %>% mutate(FID = ifelse(FID == "PFKT2312", "PFKT2312_wd", FID),
                      IID = ifelse(IID == "PFKT2312", "PFKT2312_wd", IID),
                      FID = ifelse(FID == "PFKT2334", "PFKT2334_wd", FID),
                      IID = ifelse(IID == "PFKT2334", "PFKT2334_wd", IID),
                      FID = ifelse(FID == "PFKT2377", "PFKT2377_wd", FID),
                      IID = ifelse(IID == "PFKT2377", "PFKT2377_wd", IID),
                      )

###coverage
cov <- fread("out.idepth")
cov <- cov %>% filter(MEAN_DEPTH < 10)
mis <- fread("out.imiss")
mis <- mis %>% filter(F_MISS > 0.01)

all <- all %>% mutate(tokeep = ifelse(FID %in% cov$INDV, 0, tokeep))
all <- all %>% filter(tokeep == 1)

all <- all %>% mutate(WGS = ifelse(grepl("wd", FID), "DNBSeq G400RS x30", WGS))
#filter(WGS == "DNBSeq G400RS") %>% filter(grepl("wd", FID))

write.table(all, file="ALL_sample.tsv", quote=F, col.names = T, row.names = F, sep="\t")

write.table(all$FID, file="SampleQC.tokeep.All.sample", quote=F, col.names = F, row.names = F, sep="\t")

write.table(all$FID[all$group == "Control"], file="SampleQC.tokeep.AGP3000.sample", quote=F, col.names = F, row.names = F, sep="\t")
write.table(all$FID[all$group != "Control"], file="SampleQC.tokeep.ILD.sample", quote=F, col.names = F, row.names = F, sep="\t")


cov <- fread("out.idepth") 
cov <- cov %>% inner_join(all[,c("FID", "WGS", "group")], by=c("INDV"="FID"))

cov %>% filter(group == "Control" & !grepl("FPFA", INDV)) %>% group_by(WGS) %>% 
  summarise(mean = mean(MEAN_DEPTH),
            sd = sd(MEAN_DEPTH))

cov %>% filter(!(group == "Control" & !grepl("FPFA", INDV))) %>% group_by(WGS) %>% 
  summarise(mean = mean(MEAN_DEPTH),
            sd = sd(MEAN_DEPTH))

mis <- fread("out.imiss")
mis <- mis %>% inner_join(all[,c("FID", "WGS", "group")], by=c("INDV"="FID"))

mis %>% filter(!(group == "Control" & !grepl("FPFA", INDV))) %>% group_by(WGS) %>% 
  summarise(mean = mean(F_MISS),
            sd = sd(F_MISS))

mis %>% filter(group == "Control" & !grepl("FPFA", INDV)) %>% group_by(WGS) %>% 
  summarise(mean = mean(F_MISS),
            sd = sd(F_MISS))
