setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/")

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

telomere <- fread("telomere.txt")
newbaseline <- fread("../data/IP_Diagnosis_Followup_20220916.tsv")
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

data <- data %>% left_join(newbaseline, by=c("FID"="ID"))

data <- data %>% mutate(FamilyID = case_when(IID %in% c("PFKT0751", "PFKT0901", "PFKT2176") ~ "FPFK004",
                                             IID %in% c("FPFB0001", "FPFB0002", "FPFB0005") ~ "FPFB",
                                             IID %in% c("PFKT2042", "PFKT2057", "PFKT2067") ~ "FPFK001",
                                             grepl("FPFA", IID) ~ "FPFA",
                                             IID %in% c("PFKT2486", "PFKT2531") ~ "FPFK044",
                                             TRUE ~ FID))

data <- data[!duplicated(data$FamilyID),]
data <- data %>% mutate(gene = case_when(IID %in% c("PFKT0705","PFKT0751","PFKT0901","PFKT2176","TLKU0266") ~ "MUC5B",
                                         IID %in% c("FPFB0001", "FPFB0002", "FPFB0005", 
                                                    "PFKT1358","PFKT2074","PFKT2377","PFKT2438") ~ "TERT",
                                         IID %in% c("PFKT0578") ~ "TERC",
                                         IID %in% c("PFKT2017","PFKT2042","PFKT2057","PFKT2067",
                                                    "PFKT2305","PFKT2442") ~ "SFTPA2",
                                         IID %in% c("PFKT2471") ~ "SFTPA1",
                                         IID %in% c("PFKT0568","PFKT0935","PFKT2279","TLKU9001","TLTY0055") ~ "RTEL1",
                                         IID %in% c("PFTY0190","TLTY0101") ~ "PARN",
                                         IID %in% c("TLKW0106","PFKT0160","PFKT0287","PFKT0548","PFKT0861") ~ "MUC5B",
                                         IID %in% c("TLKU0328","PFKT1222","PFKT1147","PFKT0135") ~ "TERT",
                                         IID %in% c("PFKT1679") ~ "SFTPC",
                                         IID %in% c("PFKT1808") ~ "SFTPA1",
                                         IID %in% c("PFKT0276","PFKT0427","PFKT0632","PFKT0746","PFKT0746",
                                                    "TLKU0213","PFKT1364","PFKT1290","PFKT1003") ~ "RTEL1",
                                         IID %in% c("TLKU0594","TLKU0182","TLKU0152","PFKT1679","PFKT0148") ~ "PARN",
                                         IID %in% c("TLKU0581") ~ "ABCA3",
                                         TRUE ~ "Non"))

data <- data %>% left_join(telomere, by=c("IID"="ID"))

figs <- data %>% mutate(AgeAtDiagnosis = round(as.numeric(as.Date(DateOfDiagnosis) - as.Date(DateOfBirth))/365))
figs <- figs %>% select(IID, gene, tel_sd, AgeAtDiagnosis)
figs_long <- melt(figs, measure.vars = 3:4)

library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)

png("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1BD.png", width = 700, height = 600)
figs_long %>% mutate(gene = factor(gene, levels=c("MUC5B", "RTEL1", "TERT", "PARN", "TERC","SFTPA2", "SFTPA1","SFTPC","ABCA3","Non"))) %>% 
  ggplot(aes(x=gene, y=value)) + 
  facet_wrap(~rev(variable), ncol=1, scales = "free_y") +
  geom_jitter(aes(color=gene)) +  theme_classic() +
  geom_boxplot(aes(fill=gene), color="black") + 
  theme(legend.title = element_text( size = 15),
        legend.text = element_text(size = 20, face = "italic"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 20, face = "italic", angle = 45, vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(size = 20),
        axis.title.y = element_blank()) +　
  guides(fill = FALSE, color = FALSE) +
  scale_color_manual( values = mycolors) + 
  scale_fill_manual( values = mycolors)
dev.off()
