setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/")

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
prs <- fread("PRS/IPFPRS.sscore.gz") %>% select(`#IID`, SCORE1_AVG)
colnames(prs) <- c("IID", "PRS")
data <- data %>% left_join(prs, by="IID")
data <- data %>% mutate(PRS = scale(PRS))
data <- data %>% left_join(telomere, by=c("IID"="ID"))

data <- data %>% mutate(gene = case_when(IID %in% c("FPFB0001", "FPFB0002", "FPFB0005", 
                                                    "PFKT1358","PFKT2074","PFKT2377_wd","PFKT2438") ~ "TERT",
                                         IID %in% c("PFKT0578", "NAG11988", "NAG17059", "NAG15139") ~ "TERC",
                                         IID %in% c("PFKT2017","PFKT2042","PFKT2057","PFKT2067",
                                                    "PFKT2305","PFKT2442", "7000063092", "NAG16308",
                                                    "LCAC0669", "1999934671", "NAG11096", "8000064816", "HVUM0617") ~ "SFTPA2",
                                         IID %in% c("PFKT2471", "LCAC0523", "NAG15050", "999935061", "C_AC0286", "6999935851") ~ "SFTPA1",
                                         IID %in% c("PFKT0568","PFKT0935","PFKT2279","TLKU9001","TLTY0055",
                                                    "C_AC0320", "C_AC0943", "NAG10897", "NAG13984", "NAG14098",
                                                    "NAG15170", "NAG17797", "NAG18883", "AYUM0280", "AYUM0285",
                                                    "AYUM0519", "C_AC0713", "C_AC1702", "HVUM0307", "HVUM0456",
                                                    "HVUM0733", "LCAC0130", "LCAC0848", "NAG10227", "NAG10238", 
                                                    "NAG11332", "NAG11358", "NAG12273", "NAG14396", "NAG14418",
                                                    "NAG14573", "NAG15139", "NAG15274", "NAG15375", "NAG15685",
                                                    "NAG16420", "NAG16510", "NAG16897", "NAG17298", "NAG17443",
                                                    "NAG18239", "NAG18376", "NAG18652", "NAG19135", "NAG19707",
                                                    "999934733", "999935447", "999936155", "1999936413", "2999935559",
                                                    "4000064370", "5000064002", "5000064524", "6000066332", 
                                                    "6999935203", "7000063142", "7999936051", "8000064298", 
                                                    "8999935515", "8999936427", "9000064090") ~ "RTEL1",
                                         IID %in% c("PFTY0190","TLTY0101") ~ "PARN",
                                         IID %in% c("TLKU0328","PFKT1222","PFKT1147","PFKT0135", "NAG11310") ~ "TERT",
                                         IID %in% c("PFKT1679", "C_AC0926", "LCAC0055", "LCAC0164", "NAG14089",
                                                    "C_AC0790", "C_AC1526", "NAG11645", "NAG12222", "NAG13723",
                                                    "NAG17403", "0000063154", "0000064322", "00000065928", 
                                                    "1000066082", "1999935231", "1999936699", "3000064094", 
                                                    "4000064300", "4000064390", "6000065902", "7999934285",
                                                    "7999934937", "7999936131", "7999936489", "8000063748", 
                                                    "8000064074", "9000064180") ~ "SFTPC",
                                         IID %in% c("PFKT1808") ~ "SFTPA1",
                                         IID %in% c("PFKT0276","PFKT0427","PFKT0632","PFKT0746","PFKT0746",
                                                    "TLKU0213","PFKT1364","PFKT1290","PFKT1003") ~ "RTEL1",
                                         IID %in% c("TLKU0594","TLKU0182","TLKU0152","PFKT1679","PFKT0148",
                                                    "9000064538", "NAG16860", "NAG14182","NAG13978","4999937047") ~ "PARN",
                                         IID %in% c("TLKU0581", "NAG15834") ~ "ABCA3",
                                         IID %in% c("NAG13043", "NAG14769") ~ "DKC1",
                                         IID %in% c("NAG14901", "0000064372", "1999935305") ~ "NAF1",
                                         IID %in% c("NAG14841", "NAG19685", "999933489", "4000064290", "4000065398",
                                                    "C_AC1537", "9999934031", "AYUM0400", "LCAC0130", "LCAC0207", 
                                                    "5000064640", "4999934811") ~ "STN1",
                                         IID %in% c("AYUM5117","NAG14225","HVUM0056", "5999935987", "C_AC1564", "NAG17004", "AYUM0488") ~ "TINF2",
                                         TRUE ~ "0Non"))

data <- data %>% mutate(telomere_gene = case_when(gene %in% c("TERT", "RTEL1", "PARN", "TERC", "NAF1", "TINF2", "DKC1", "STN1") ~ 1,
                                                  TRUE ~ 0))

data <- data %>% mutate(PRSquantile = ntile(PRS, 5))
data <- data %>% mutate(TSquantile = ntile(tel_sd, 5))
data <- data %>% mutate(smoking = case_when(smoking == "never" ~ "0never",
                                            smoking == "ever" ~ "ever"))
data <- data %>% arrange(desc(Status))

toremove <- fread("kingunrelated_toberemoved.txt", header = F)
toremove <- toremove %>% filter(!(grepl("PF|TL", V2)))
data1 <- data %>% filter(!(IID %in% toremove$V2))


data1 <- data1 %>% filter(PCAoutlier == FALSE)

fpf <- data1 %>% filter(Status %in% c("FPF (IPF)", "FPF (non-IPF)")) 
fpf <- fpf %>% mutate(diagnosis = case_when(Status == "FPF (IPF)" ~ "IPF",
                                            IID %in% c("PFKT0109") ~ "ARS-ILD",
                                            IID %in% c("PFKT2067") ~ "CTD-ILD",
                                            IID %in% c("PFKT0578", "PFKT0568", "PFKT0262", "PFKT2222") ~ "HP",
                                            IID %in% c("PFKT0578", "PFKT0568", "PFKT0262", "PFKT2222") ~ "HP",
                                            IID %in% c("PFKT2176") ~ "MPA-ILD",
                                            IID %in% c("PFKT1304") ~ "NSIP",
                                            IID %in% c("PFKT1358", "PFKT2074", "PFKT2438", "PFKT2279", 
                                                       "PFKT2282", "PFKT2405") ~ "PPFE",
                                            IID %in% c("PFKT2377_wd", "PFKT2513", "PFKT2546", "TLKU0699") ~ "SjS-ILD",
                                            IID %in% c("PFKT2293") ~ "SSc-ILD",
                                            IID %in% c("PFKT2017", "PFKT2305", "PFKT2442", "PFKT0112", 
                                                       "PFKT2312_wd", "PFKT2334_wd", "TLKU0535", "PFKT2492", 
                                                       "PFKT2500", "PFKT2525", "PFKT2531") ~ "Unclassifiable")
                                            
)

data1 <- data1 %>% mutate(FPF = case_when(Status == "FPF (IPF)" ~ 1,
                                          Status == "Control" ~ 0))
fpf <- data1 %>% drop_na(PRS, telomere_gene, tel_sd, age.x, geneticSex.x, smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
fpf <- fpf[!duplicated(fpf$FamilyID),]

LM <- glm(FPF ~ PRS + telomere_gene + tel_sd + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, family="binomial", data=fpf)


out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(variable = rownames(out))
out <- out %>% mutate(Estimate = case_when(variable == "tel_sd" ~ -1*Estimate,
                                           variable == "age.x" ~ 10*Estimate,
                                           TRUE  ~ Estimate))
out <- out %>% mutate(OR = exp(Estimate),
                      LL = exp(Estimate + qnorm(0.025)*Std..Error),
                      UL = exp(Estimate + qnorm(0.975)*Std..Error), 
                      p=Pr...z..,
                      Ncases = table(fpf$FPF)[2],
                      Ncontrols = table(fpf$FPF)[1])

out <- out %>% mutate(variable = case_when(variable %in% c("PRS", "prs") ~ "1SD increase in PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("tel_sd", "TS_sd")~ "1SD decrease in age/sex adjusted LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))
out <- out %>% drop_na(variable)
out <- out %>% select(variable, OR, LL, UL, p, Ncases, Ncontrols)

out %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/FPFIPF.susceptibility.assoc.xlsx")


data1 <- data1 %>% mutate(HighPRS = ifelse(PRSquantile == 5, 1, 0),
                      ShortTS = ifelse(TSquantile == 1, 1, 0))
#data1$ShortTS[is.na(data1$ShortTS)] <- 0
fpf <- data1 %>% drop_na(PRS, telomere_gene, tel_sd, age.x, geneticSex.x, smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
fpf <- fpf[!duplicated(fpf$FamilyID),]

LM <- glm(FPF ~ HighPRS + telomere_gene + ShortTS + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, family="binomial", data=fpf)

summary(LM)
out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(variable = rownames(out))
out <- out %>% mutate(Estimate = case_when(variable == "tel_sd" ~ -1*Estimate,
                                           variable == "age.x" ~ 10*Estimate,
                                           TRUE  ~ Estimate))
out <- out %>% mutate(OR = exp(Estimate),
                      LL = exp(Estimate + qnorm(0.025)*Std..Error),
                      UL = exp(Estimate + qnorm(0.975)*Std..Error), 
                      p=Pr...z..,
                      Ncases = table(fpf$FPF)[2],
                      Ncontrols = table(fpf$FPF)[1])

out <- out %>% mutate(variable = case_when(variable %in% c("HighPRS") ~ "High IPF-PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("ShortTS")~ "Short LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))
out <- out %>% drop_na(variable)
out_IPF <- out %>% mutate(outcome = "familial IPF")
out <- out %>% select(variable, OR, LL, UL, p, Ncases, Ncontrols)

out %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/categorized_FPFIPF.susceptibility.assoc.xlsx")

FPF <- fpf %>% filter(FPF == 1)
FPF <- FPF %>% mutate(Risk = case_when(telomere_gene == 1 & HighPRS == 1 & ShortTS == 1 ~ "Rare telomere gene carrier AND high PRS AND short LTL",
                                       telomere_gene == 1 & HighPRS == 0 & ShortTS == 1 ~ "Rare telomere gene carrier AND low PRS AND short LTL",
                                       telomere_gene == 1 & HighPRS == 1 & ShortTS == 0 ~ "Rare telomere gene carrier AND high PRS AND long LTL",
                                       telomere_gene == 1 & HighPRS == 0 & ShortTS == 0 ~ "Rare telomere gene carrier AND low PRS AND long LTL",
                                       telomere_gene == 0 & HighPRS == 1 & ShortTS == 1 ~ "Rare telomere gene non-carrier AND high PRS AND short LTL",
                                       telomere_gene == 0 & HighPRS == 0 & ShortTS == 1 ~ "Rare telomere gene non-carrier AND low PRS AND short LTL",
                                       telomere_gene == 0 & HighPRS == 1 & ShortTS == 0 ~ "Rare telomere gene non-carrier AND high PRS AND long LTL",
                                       telomere_gene == 0 & HighPRS == 0 & ShortTS == 0 ~ "Rare telomere gene non-carrier AND low PRS AND long LTL",
))

table(FPF$Risk)

data1 <- data1 %>% mutate(FPF = case_when(Status == "FPF (non-IPF)" ~ 1,
                                          Status == "Control" ~ 0))
fpf <- data1 %>% drop_na(PRS, telomere_gene, tel_sd, age.x, geneticSex.x, smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
fpf <- fpf[!duplicated(fpf$FamilyID),]

LM <- glm(FPF ~ PRS + telomere_gene + tel_sd + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, family="binomial", data=fpf)

FPF <- fpf %>% filter(FPF == 1)
FPF <- FPF %>% mutate(Risk = case_when(telomere_gene == 1 & HighPRS == 1 & ShortTS == 1 ~ "Rare telomere gene carrier AND high PRS AND short LTL",
                                       telomere_gene == 1 & HighPRS == 0 & ShortTS == 1 ~ "Rare telomere gene carrier AND low PRS AND short LTL",
                                       telomere_gene == 1 & HighPRS == 1 & ShortTS == 0 ~ "Rare telomere gene carrier AND high PRS AND long LTL",
                                       telomere_gene == 1 & HighPRS == 0 & ShortTS == 0 ~ "Rare telomere gene carrier AND low PRS AND long LTL",
                                       telomere_gene == 0 & HighPRS == 1 & ShortTS == 1 ~ "Rare telomere gene non-carrier AND high PRS AND short LTL",
                                       telomere_gene == 0 & HighPRS == 0 & ShortTS == 1 ~ "Rare telomere gene non-carrier AND low PRS AND short LTL",
                                       telomere_gene == 0 & HighPRS == 1 & ShortTS == 0 ~ "Rare telomere gene non-carrier AND high PRS AND long LTL",
                                       telomere_gene == 0 & HighPRS == 0 & ShortTS == 0 ~ "Rare telomere gene non-carrier AND low PRS AND long LTL",
))

table(FPF$Risk)

out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(variable = rownames(out))
out <- out %>% mutate(Estimate = case_when(variable == "tel_sd" ~ -1*Estimate,
                                           variable == "age.x" ~ 10*Estimate,
                                           TRUE  ~ Estimate))
out <- out %>% mutate(OR = exp(Estimate),
                      LL = exp(Estimate + qnorm(0.025)*Std..Error),
                      UL = exp(Estimate + qnorm(0.975)*Std..Error), 
                      p=Pr...z..,
                      Ncases = table(fpf$FPF)[2],
                      Ncontrols = table(fpf$FPF)[1])

out <- out %>% mutate(variable = case_when(variable %in% c("PRS", "prs") ~ "1SD increase in PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("tel_sd", "TS_sd")~ "1SD decrease in age/sex adjusted LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))
out <- out %>% drop_na(variable)
out <- out %>% select(variable, OR, LL, UL, p, Ncases, Ncontrols)

out %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/FPFnonIPF.susceptibility.assoc.xlsx")


data1 <- data1 %>% mutate(HighPRS = ifelse(PRSquantile == 5, 1, 0),
                          ShortTS = ifelse(TSquantile == 1, 1, 0))

LM <- glm(FPF ~ HighPRS + telomere_gene + ShortTS + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, family="binomial", data=fpf)

summary(LM)
out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(variable = rownames(out))
out <- out %>% mutate(Estimate = case_when(variable == "tel_sd" ~ -1*Estimate,
                                           variable == "age.x" ~ 10*Estimate,
                                           TRUE  ~ Estimate))
out <- out %>% mutate(OR = exp(Estimate),
                      LL = exp(Estimate + qnorm(0.025)*Std..Error),
                      UL = exp(Estimate + qnorm(0.975)*Std..Error), 
                      p=Pr...z..,
                      Ncases = table(fpf$FPF)[2],
                      Ncontrols = table(fpf$FPF)[1])

out <- out %>% mutate(variable = case_when(variable %in% c("HighPRS") ~ "High IPF-PRS",
                                           variable %in% c("telomere_gene", "tel_gene") ~ "Rare telomere variants carrier",
                                           variable %in% c("ShortTS")~ "Short LTL",
                                           variable %in% c("age.x", "AGE") ~ "10 year increase in age at recruitment",
                                           variable %in% c("geneticSex.xMale", "SEX.y") ~ "Male",
                                           variable %in% c("smokingever", "smokingEver") ~ "Ever smoker"))
out <- out %>% drop_na(variable)
out_nonIPF <- out %>% mutate(outcome = "nonIPF FPF")

out <- out %>% select(variable, OR, LL, UL, p, Ncases, Ncontrols)

out %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/categorized_FPFnonIPF.susceptibility.assoc.xlsx")


IPF <- readxl::read_excel("../repo/IPF_GWAS/Tables/categorized_FPFIPF.susceptibility.assoc.xlsx") %>% mutate(outcome = "Familial IPF")
nonIPF <- readxl::read_excel("../repo/IPF_GWAS/Tables/categorized_FPFnonIPF.susceptibility.assoc.xlsx") %>% mutate(outcome = "non-IPF FPF")

ALL <- bind_rows(IPF, nonIPF)
library(forcats)
ALL %>% mutate(name = fct_relevel(variable, "Rare telomere variants carrier", "High IPF-PRS", "Short LTL",
                                  "10 year increase in age at recruitment", "Male", "Ever smoker")) %>%
  mutate(outcome1 = fct_relevel(outcome, "non-IPF FPF", "Familial IPF")) %>% 
  ggplot(aes(x=outcome1, y=OR, ymin=LL, ymax=UL, color=outcome1)) +
  geom_pointrange(aes(col=outcome1), lwd=0.8) + geom_hline(aes(fill=outcome1), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + geom_text(aes(label=round(OR,1), y=OR, col=outcome1, hjust = 0.5, vjust = 1.5), size=5) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=outcome1), width=0.1, cex=1) + #ylim(0.8,3.5)+
  facet_wrap(~name,  strip.position = 'left', nrow = 10) + theme_minimal() +scale_y_log10() +
  scale_color_manual(labels = c("non-IPF FPF", "Familial IPF"), values=c("#fc8d59", "#91bfdb")) + labs(color="Cohort") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))

ggsave(paste0("../repo/IPF_GWAS/Figures/Fig5a.png"), width=25, height=15, units = "cm", dpi=300)



FPF <- data1 %>% filter(FPF == 1)
FPF <- FPF %>% mutate(Risk = case_when(telomere_gene == 1 & HighPRS == 1 & ShortTS == 1 ~ "Rare telomere gene carrier AND high PRS AND short LTL",
                                       telomere_gene == 1 & HighPRS == 0 & ShortTS == 1 ~ "Rare telomere gene carrier AND low PRS AND short LTL",
                                       telomere_gene == 1 & HighPRS == 1 & ShortTS == 0 ~ "Rare telomere gene carrier AND high PRS AND long LTL",
                                       telomere_gene == 1 & HighPRS == 0 & ShortTS == 0 ~ "Rare telomere gene carrier AND low PRS AND long LTL",
                                       telomere_gene == 0 & HighPRS == 1 & ShortTS == 1 ~ "Rare telomere gene non-carrier AND high PRS AND short LTL",
                                       telomere_gene == 0 & HighPRS == 0 & ShortTS == 1 ~ "Rare telomere gene non-carrier AND low PRS AND short LTL",
                                       telomere_gene == 0 & HighPRS == 1 & ShortTS == 0 ~ "Rare telomere gene non-carrier AND high PRS AND long LTL",
                                       telomere_gene == 0 & HighPRS == 0 & ShortTS == 0 ~ "Rare telomere gene non-carrier AND low PRS AND long LTL",
))

table(FPF$Risk)

FPF %>% filter(HighPRS == 1)



ALL <- bind_rows(out_IPF, out_nonIPF)
cov <- unique(ALL$variable)
for(j in c(1:length(cov))){
  tmp <- ALL %>% filter(variable == cov[j]) 
  m1 <- meta::metagen(Estimate,
                      Std..Error,
                      data=tmp,
                      studlab=paste(outcome),
                      comb.fixed = TRUE,
                      comb.random = TRUE,
                      prediction = FALSE,
                      sm="OR")
  ALL$hetero_p[ALL$variable == cov[j]] <- m1$pval.Q
}

ALL <- ALL %>% select(variable, OR, LL, UL, p, hetero_p, Ncases, Ncontrols, outcome)
ALL %>% openxlsx::write.xlsx("../repo/IPF_GWAS/Tables/categorized_FPF.susceptibility.assoc.xlsx")



data2 <- data %>% filter(PCAoutlier == FALSE)
data2 <- data2 %>% filter(Status %in% c("FPF (IPF)", "FPF (non-IPF)") | grepl("FPFA", IID))
data2 <- data2 %>% drop_na(FID)
data2 %>% filter(is.na(telomere_gene))
data2 %>% filter(gene != "0Non") %>% select(IID, gene) %>% View()

##survival
# clinical <- read.csv("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/clinical/clinical20230605.csv")
# clinical1 <- clinical %>% select(`難プラID`, `生存確認コード`, `最終生存確認日`)
# clinical1 <- clinical1 %>% mutate(DateOfDeath = as.Date(`最終生存確認日`, format="%Y/%m/%d"))
# clinical1 <- clinical1 %>% drop_na(DateOfDeath)#
# clinical1 <- clinical1 %>% group_by(`難プラID`) %>%
#   filter(DateOfDeath == max(DateOfDeath)) %>% ungroup() %>%
#   distinct(`難プラID`, .keep_all=TRUE)
# clinical1 <- clinical1 %>% select(-`最終生存確認日`)
# clinical2 <- clinical %>%  select(`難プラID`, `生存確認コード`, `死亡日`)
# clinical2 <- clinical2 %>% mutate(DateOfDeath = as.Date(`死亡日`, format="%Y/%m/%d"))
# clinical2 <- clinical2 %>% drop_na(DateOfDeath)#
# clinical2 <- clinical2 %>% group_by(`難プラID`) %>%
#   filter(DateOfDeath == max(DateOfDeath)) %>% ungroup() %>%
#   distinct(`難プラID`, .keep_all=TRUE)
# clinical2 <- clinical2 %>% select(-`死亡日`)
# 
# clinical <- bind_rows(clinical1, clinical2)
# clinical <- clinical %>% arrange(desc(生存確認コード))
# clinical <- clinical[!duplicated(clinical$難プラID),]
# 
# clinical <- clinical %>% mutate(death = ifelse(生存確認コード %in% c(1, 999), 0, 1))
# mapping <- read.csv("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/clinical/mapping20230605.csv")
# mapping <- mapping %>% select(`難プラID`, `間質性肺炎の個別ID`)
# clinical <- clinical %>% left_join(mapping, by="難プラID")
# #/home/pub/data/projects/Interstitial_Pneumonia/Sample/telomere/IP_telomere検体一覧_20220304.xlsx
# #/home/pub/data/projects/Interstitial_Pneumonia/Sample/IP検体一覧20230509.xlsx
# tmp1 <- readxl::read_excel("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/data/IP検体一覧20230509.xlsx")
# tmp1 <- tmp1 %>% select(`検体 ID`, `CGM ID`)
# tmp1 <- unique(tmp1)
# colnames(tmp1) <- c("ID1", "ID2")
# tmp2 <- readxl::read_excel("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/data/IP_telomere検体一覧_20220304.xlsx")
# tmp2 <- tmp2 %>% select(`症例ID`, `CGM patient ID`)
# tmp2 <- unique(tmp2)
# colnames(tmp2) <- c("ID1", "ID2")
# tmp <- bind_rows(tmp1, tmp2)
# tmp <- unique(tmp)
# tmp %>% filter(nchar(ID1) != 7) %>% View()
# tmp <- tmp %>% mutate(ID1 = case_when(nchar(ID1) == 6 ~ paste0(substr(ID1, 1,3),"0",substr(ID1, 4,6)),
#                                       grepl("-", ID1) ~ substr(ID1, 1, 7),
#                                       TRUE ~ ID1))
# tmp <- tmp %>% mutate(ID2 = case_when(grepl("_", ID2) ~ substr(ID2, 1, 8),
#                                       TRUE ~ ID2))
# 
# clinical <- clinical %>% mutate(間質性肺炎の個別ID = ifelse(間質性肺炎の個別ID == "PFK225", "PFK0225", 間質性肺炎の個別ID))
# clinical <- clinical %>% left_join(tmp, by=c("間質性肺炎の個別ID"="ID1"))
# 
# diagnosis <- read.csv("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/clinical/diagnosis20230613.csv")
# diagnosis <- diagnosis %>% select(`難プラID`, 難病診断年月.最初の診断.)
# diagnosis <- diagnosis %>% mutate(難病診断年月.最初の診断. = paste0(難病診断年月.最初の診断.,"/01"))
# diagnosis <- diagnosis %>% mutate(DateOfDiagnosis = as.Date(難病診断年月.最初の診断., format="%Y/%m/%d"))
# diagnosis <- diagnosis %>% select(-難病診断年月.最初の診断.) %>% unique()
# clinical <- clinical %>% left_join(diagnosis, by="難プラID")
# clinical <- clinical %>% mutate(ID2 = case_when(間質性肺炎の個別ID == "PFK0071" ~ "PFKT0710",
#                                                 間質性肺炎の個別ID == "PFK0077" ~ "PFKT0773",
#                                                 間質性肺炎の個別ID == "PFK0116" ~ "TLKU9018",
#                                                 間質性肺炎の個別ID == "PFK0150" ~ "PFKT1506",
#                                                 間質性肺炎の個別ID == "PFK0168" ~ "PFKK1684",
#                                                 間質性肺炎の個別ID == "PFK0129" ~ "PFKT1290",
#                                                 TRUE ~ ID2))
# head(clinical)
# clinical <- clinical %>% select(ID2, DateOfDiagnosis, death, DateOfDeath)
# clinical <- clinical %>% drop_na()
# clinical <- clinical %>% mutate(Start = as.Date(DateOfDiagnosis),
#                                 End_death = as.Date(DateOfDeath),
#                                 ID = ID2) %>%
#   dplyr::select(ID, death, Start, End_death)

follow <- openxlsx::read.xlsx("IP_Diagnosis_Followup20221018_TN.xlsx")
follow <- follow %>% dplyr::select(ID, DateOfDiagnosis, death, DateOfDeath, death_date)

follow <- follow %>% mutate(Start = as.Date(DateOfDiagnosis),
                            End_death = case_when(!is.na(death_date) ~ as.Date(as.numeric(death_date), origin = "1899-12-30"),
                                                  TRUE ~ as.Date(DateOfDeath)))

follow <- follow %>% dplyr::select(ID, death, Start, End_death)
follow <- follow %>% drop_na()

follow <- bind_rows(clinical, follow)
follow <- follow %>% arrange(desc(death))
follow <- follow[!duplicated(follow$ID),]

library(survival)
data1 <- data1 %>% mutate(HighPRS = ifelse(PRSquantile == 5, 1, 0),
                          ShortTS = ifelse(TSquantile == 1, 1, 0))
ipf <- data1 %>% filter(IPF == 1)
ipf <- ipf[!duplicated(ipf$FamilyID),]

follow1 <- follow %>% right_join(ipf, by=c("ID" = "IID"))
follow1 <- follow1 %>% mutate(time = (End_death - Start)/30)
follow1 <- follow1 %>% drop_na(time)

follow1 <- follow1 %>% mutate(death = ifelse(death.x == 2, NA, death.x))
follow1 <- follow1 %>% select(time, death, PRS, telomere_gene, tel_sd, age.x, geneticSex.x,
                              smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
follow1 <- follow1[complete.cases(follow1),]
quantile(follow1$time)
ip.cox <- coxph(Surv(time, death) ~ PRS + telomere_gene + tel_sd + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, data = follow1)
summary(ip.cox)

survfit(Surv(time, death) ~ 1, data = follow1)

out1 <- data.frame(summary(ip.cox)$coefficients)
out1 <- out1 %>% dplyr::select(Pr...z.., coef, se.coef.)
out2 <- data.frame(summary(ip.cox)$conf.int)
out2 <- out2 %>% dplyr::select(-exp..coef.)

out <- bind_cols(out1, out2)
out <- out %>% mutate(variable = rownames(out))
out <- out %>% mutate(Ncases = 34, Ncontrols = 89)

out %>% write.table("IPF.survival.assoc.tsv", sep="\t", quote=F, col.names = T, row.names = F)

follow1 <- follow %>% right_join(ipf, by=c("ID" = "IID"))
follow1 <- follow1 %>% mutate(time = (End_death - Start)/30)
follow1 <- follow1 %>% drop_na(time)

follow1 <- follow1 %>% mutate(death = ifelse(death.x == 2, NA, death.x))
follow1 <- follow1 %>% select(time, death, HighPRS, telomere_gene, ShortTS, age.x, geneticSex.x,
                              smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
follow1 <- follow1[complete.cases(follow1),]

ip.cox <- coxph(Surv(time, death) ~ HighPRS + telomere_gene + ShortTS + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, data = follow1)
summary(ip.cox)

survfit(Surv(time, death) ~ telomere_gene, data = follow1)

out1 <- data.frame(summary(ip.cox)$coefficients)
out1 <- out1 %>% dplyr::select(Pr...z.., coef, se.coef.)
out2 <- data.frame(summary(ip.cox)$conf.int)
out2 <- out2 %>% dplyr::select(-exp..coef.)

out <- bind_cols(out1, out2)
out <- out %>% mutate(variable = rownames(out))
out <- out %>% mutate(Ncases = 34, Ncontrols = 89)

out %>% write.table("category_IPF.survival.assoc.tsv", sep="\t", quote=F, col.names = T, row.names = F)

