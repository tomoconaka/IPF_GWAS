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

data <- data %>% inner_join(telomere, by=c("IID"="ID"))

data <- data %>% mutate(FamilyID = case_when(IID %in% c("PFKT0751", "PFKT0901", "PFKT2176") ~ "FPFK004",
                                             IID %in% c("FPFB0001", "FPFB0002", "FPFB0005") ~ "FPFB",
                                             IID %in% c("PFKT2042", "PFKT2057", "PFKT2067") ~ "FPFK001",
                                             grepl("FPFA", IID) ~ "FPFA",
                                             IID %in% c("PFKT2486", "PFKT2531") ~ "FPFK044",
                                             TRUE ~ FID))

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


fig3 <- data %>% filter(DNB == FALSE)
fig3 <- fig3 %>% mutate(status = case_when(Status %in% c("FPF (IPF)", "sporadic IPF") ~ "IPF", 
                                           Status == "Control" ~ "0Control",
                                           TRUE ~ "Other_ILD"))

LM <- glm(tel_sd ~ status, dat=fig3, family="gaussian")
summary(LM)
out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(cohort = "Kyoto",
                      LL = Estimate - qnorm(0.975)*Std..Error,
                      UL = Estimate + qnorm(0.975)*Std..Error,
                      Status = gsub("status", "", rownames(out)))
out$N <- as.numeric(table(fig3$status)[-1])
Kyoto <- out %>% select(Status, Estimate, LL, UL, Pr...t.., cohort, N)
kyoto_add <- c("Control", 0,0,0,1, "Kyoto",3028)
kyoto_add <- data.frame(t(kyoto_add))
colnames(kyoto_add) <- colnames(Kyoto)
kyoto_add[,c(2:5,7)] <- as.numeric(kyoto_add[,c(2:5,7)])
Kyoto <- bind_rows(Kyoto, kyoto_add)

UKB <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/candidate/TS_UKB_summary.tsv")
ukb_add <- c("Control", 0,0,0,1, "UKB",423847)
ukb_add <- data.frame(t(ukb_add))
colnames(ukb_add) <- colnames(UKB)
ukb_add[,c(2:5,7)] <- as.numeric(ukb_add[,c(2:5,7)])
UKB <- bind_rows(UKB, ukb_add)
out <- bind_rows(Kyoto, UKB)

library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(11)
library(gridExtra)
png("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1_Kyoto.png", width = 700, height = 200)
p1 <- Kyoto %>% mutate(Status = factor(Status, levels=c("Other_ILD", "IPF","Control"))) %>% 
  ggplot(aes(x=Status, y=Estimate, ymin=LL, ymax=UL, color=Status)) +
  geom_pointrange(aes(col=Status), lwd=0.8) + geom_hline(aes(fill=Status), yintercept =0, linetype=2) +
  xlab("") + ylab("") +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=Status), width=0.1, cex=1) + ylim(min(out$LL),max(out$UL))+
  theme_minimal() +
  scale_color_manual(values = mycolors[c(9,1,11)]) + labs(color="Kyoto Status") +
  scale_x_discrete(
    labels = c("Control"=paste0("Control [ref]\n(N=",Kyoto$N[Kyoto$Status == "Control"],")"),
               "Other_ILD"=paste0("Non-IPF ILD\n(N=",Kyoto$N[Kyoto$Status == "Other_ILD"],")"),
               "IPF"=paste0("IPF\n(N=",Kyoto$N[Kyoto$Status == "IPF"],")"))) +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_text(size=15,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_blank(),
        legend.position="none",
        legend.title = element_blank())+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))
p1
dev.off()

png("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1_UKB.png", width = 700, height = 500)

p2 <- UKB %>% mutate(Status = factor(Status, levels=c("ILD_others","SLE_ILD","DM_ILD","SSc_ILD","ANCA_ILD","SjS_ILD","RA_ILD", "IPF","Control"))) %>% 
  ggplot(aes(x=Status, y=Estimate, ymin=LL, ymax=UL, color=Status)) +
  geom_pointrange(aes(col=Status), lwd=0.8) + geom_hline(aes(fill=Status), yintercept =0, linetype=2) +
  xlab("") + ylab("Standardised telomere length difference (95% Confidence Interval)") +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=Status), width=0.1, cex=1) + ylim(min(out$LL),max(out$UL))+
  theme_minimal() +
  scale_color_manual(values = rev(mycolors[c(11,1:8)])) + labs(color="Status") +
  scale_x_discrete(
    labels = c("Control"=paste0("Control [ref]\n(N=",UKB$N[UKB$Status == "Control"],")"),
               "ANCA_ILD"=paste0("ANCA related ILD\n(N=",UKB$N[UKB$Status == "ANCA_ILD"],")"),
               "DM_ILD"=paste0("DM related ILD\n(N=",UKB$N[UKB$Status == "DM_ILD"],")"),
               "ILD_others"=paste0("Other ILD\n(N=",UKB$N[UKB$Status == "ILD_others"],")"),
               "IPF"=paste0("IPF\n(N=",UKB$N[UKB$Status == "IPF"],")"),
               "RA_ILD"=paste0("RA related ILD\n(N=",UKB$N[UKB$Status == "RA_ILD"],")"),
               "SLE_ILD"=paste0("SLE related ILD\n(N=",UKB$N[UKB$Status == "SLE_ILD"],")"),
               "SSc_ILD"=paste0("SSc related ILD\n(N=",UKB$N[UKB$Status == "SSc_ILD"],")"),
               "SjS_ILD"=paste0("SjS related ILD\n(N=",UKB$N[UKB$Status == "SjS_ILD"],")"))) +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_text(size=15,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_blank(), 
        legend.title = element_blank(), 
        legend.position="none",
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))
p2
dev.off()
