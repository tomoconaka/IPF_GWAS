setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/")

library(data.table)
library(tidyr)
library(dplyr)

data <- fread("ALL_sample_PCs.tsv")

sex <- fread("chr.idxstats.count.txt")
sex <- sex %>% mutate(Xratio = chrX/mcount,
                      Yratio = chrY/mcount)
sex <- sex %>% mutate(geneticSex = ifelse(Xratio < 0.0275, "Male", "Female"))

agp3000 <- fread("AGP3K_agesexsmoking20210819.txt")
agp3000 <- agp3000 %>% select(ID, age, `sex(2=female)`, `smoking (2 categories)`)
colnames(agp3000) <- c("ID", "age", "sex", "smoking")
ipf <- fread("../data/IPsample20220601.list")
ipf <- ipf %>% select(`CGM ID (二次)`, age, sex, SmokingStatus)
ipf <- ipf %>% mutate(sex = ifelse(sex == 1, 2, 1), 
                      smoking = case_when(SmokingStatus %in% c("ex", "cur", "current", "current?","ex?")  ~ "ever",
                                          SmokingStatus == "never" ~ "never"))
ipf <- ipf %>% select(-SmokingStatus)
colnames(ipf) <- c("ID", "age", "sex", "smoking")
ipf <- ipf %>% mutate(ID = case_when(ID == "PFKT2312" ~ "PFKT2312_wd",
                                     ID == "PFKT2334" ~ "PFKT2334_wd",
                                     ID == "PFKT2377" ~ "PFKT2377_wd",
                                     TRUE ~ ID))
ipf$age <- as.numeric(ipf$age)

sex <- sex %>% filter(id %in% data$FID)

sex <- sex %>% mutate(geneticSex = ifelse(Xratio < 0.04, "Male", "Female"))

sex1 <- sex %>% select(id, geneticSex)
data <- data %>% inner_join(sex1, by=c("FID"="id"))


dat <- bind_rows(agp3000, ipf)

data <- data %>% left_join(dat, by=c("FID"="ID"))

data %>% filter(geneticSex == "Male" & sex.x == "Female")
data %>% filter(geneticSex == "Female" & sex.x == "Male")

data <- data %>% mutate(IPF = case_when(group == "sIPF" ~ 1,
                                        group == "fIPF" ~ 1,
                                        group == "Control" ~ 0),
                        ILD = case_when(group %in% c("fIPF", "FPF", "sIPF")  ~1,
                                        TRUE ~ 0),
                        FPF = case_when(group %in% c("fIPF", "FPF")  ~1,
                                        group == "Control" ~ 0),
                        sIPF = case_when(group %in% c("sIPF")  ~1,
                                         group == "Control" ~ 0))


data1 <- data %>% group_by(group) %>% 
  mutate(age = ifelse(is.na(age), mean(age, na.rm=T), age))


pheno <- data1 %>% filter(PCAoutlier == FALSE)
pheno <- pheno %>% mutate(platform = ifelse(grepl("DNB", WGS), "MGI", "Illumina"))
pheno <- pheno %>% select(FID, IID, IPF, ILD, FPF, sIPF, age, geneticSex, PC1:PC10, platform)

pheno <- pheno[,-1]
write.table(pheno, "pheno.txt", quote=F, col.names = T, row.names = F, sep="\t")
table(pheno$sIPF)


data %>% filter(geneticSex == "Male") %>% dim()
