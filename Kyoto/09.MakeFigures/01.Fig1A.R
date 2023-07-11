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

data <- data %>% inner_join(newbaseline, by="ID")

fig1 <- data %>% filter(Status != "Control")
fig1 <- fig1 %>% mutate(FamilyID = case_when(IID %in% c("PFKT0751", "PFKT0901", "PFKT2176") ~ "FPFK004",
                                     IID %in% c("FPFB0001", "FPFB0002", "FPFB0005") ~ "FPFB",
                                     IID %in% c("PFKT2042", "PFKT2057", "PFKT2067") ~ "FPFK001",
                                     grepl("FPFA", IID) ~ "FPFA",
                                     IID %in% c("PFKT2486", "PFKT2531") ~ "FPFK044",
                                     TRUE ~ FID))

fig1 <- fig1[!duplicated(fig1$FamilyID),]
fig1 <- fig1 %>% mutate(gene = case_when(IID %in% c("PFKT0705","PFKT0751","PFKT0901","PFKT2176","TLKU0266") ~ "MUC5B",
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

case1 <- fig1 %>% filter(Status != "sporadic IPF")
tmp1 <- binom.test(sum(case1$gene != "Non"), dim(case1)[1])
case1_gene <- case1 %>% 
  group_by(gene) %>% 
  summarise(Frequency = sum(!is.na(IID))/dim(case1)[1]) %>% ungroup() %>%
  distinct(gene, .keep_all=TRUE) %>% 
  mutate(group = paste0("FPF\n",sum(case1$gene != "Non"),"/",dim(case1)[1],
                        "\n",round(sum(case1$gene != "Non")/dim(case1)[1]*100, 1),"[",round(tmp1$conf.int[1]*100, 1),"; ",round(tmp1$conf.int[2]*100, 1),"] %"))
case2 <- fig1 %>% filter(Status == "sporadic IPF")
tmp2 <- binom.test(sum(case2$gene != "Non"), dim(case2)[1])
case2_gene <- case2 %>% 
  group_by(gene) %>% 
  summarise(Frequency = sum(!is.na(IID))/dim(case2)[1]) %>% ungroup() %>%
  distinct(gene, .keep_all=TRUE) %>% 
  mutate(group = paste0("Sporadic IPF\n",sum(case2$gene != "Non"),"/",dim(case2)[1],
                        "\n",round(sum(case2$gene != "Non")/dim(case2)[1]*100, 1)," [",round(tmp2$conf.int[1]*100, 1),"; ",round(tmp2$conf.int[2]*100, 1),"] %"))


case <- bind_rows(case1_gene, case2_gene)
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(10)

png("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1A.png", width = 700, height = 300)
case %>% mutate(gene = factor(gene, levels=c("MUC5B", "RTEL1", "TERT", "PARN", "TERC","SFTPA2", "SFTPA1","SFTPC","ABCA3","Non"))) %>% 
  ggplot(aes(x="", y=Frequency, fill=gene)) +
  geom_bar(stat="identity", width=1, color="white",position = position_stack(reverse = TRUE)) +
  coord_polar("y", start=0) + facet_wrap(~group) +
  scale_fill_manual("Individuals\nwith putative pathogenic variant\nin candidate genes", values = mycolors) + 
  theme_void() + theme(legend.title = element_text( size = 15),
                       legend.text = element_text(size = 20, face = "italic"),
                       axis.title.x = element_blank(),
                       axis.text.x  = element_blank(),
                       axis.title.y = element_blank(),
                       strip.text.x = element_text(size = 20))
dev.off()





ipf <- ipf %>% mutate(DateOfBirth = as.Date(birthday, format="%d/%m/%Y"),
                      DateOfRecruitment = as.Date(recruteday, format="%d/%m/%Y"))
ipf <- ipf %>% select(`CGM ID (二次)`, DateOfBirth, DateOfRecruitment)
colnames(ipf)[1] <- c("ID")

clinical <- fread("/home/tnakanishi/gpfs1/01.IPFWGS/f2/01.basicQC/01.SampleQC/07.clinical/IPbaseline.txt")

data <- bind_rows(ipf, clinical)



