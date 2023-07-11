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
data <- data %>% inner_join(telomere, by=c("IID"="ID"))

data <- data %>% mutate(gene = case_when(IID %in% c("FPFB0001", "FPFB0002", "FPFB0005", 
                                                    "PFKT1358","PFKT2074","PFKT2377","PFKT2438") ~ "TERT",
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
data1 <- data[!duplicated(data$FamilyID),]
toremove <- toremove %>% filter(!(grepl("PF|TL", V2)))
data1 <- data1 %>% filter(!(IID %in% toremove$V2))
data1 <- data1 %>% filter(PCAoutlier == FALSE)
data1 <- data1 %>% mutate(IPF = case_when(grepl("IPF", Status) ~ 1,
                                          Status == "Control" ~ 0))

tmp <- data1 %>% select(IPF, PRS, telomere_gene, tel_sd, age.x, 
                        geneticSex.x, smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
tmp <- tmp[complete.cases(tmp),] 
table(tmp$IPF)
LM <- glm(IPF ~ PRS + telomere_gene + tel_sd + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, family="binomial", data=data1)

out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(variable = rownames(out))

out %>% write.table("IPF.susceptibility.assoc.tsv", sep="\t", quote=F, col.names = T, row.names = F)

data1 <- data1 %>% mutate(HighPRS = ifelse(PRSquantile == 5, 1, 0),
                          ShortTS = ifelse(TSquantile == 1, 1, 0))

tmp <- data1 %>% select(IPF, HighPRS, telomere_gene, ShortTS, age.x, 
                        geneticSex.x, smoking, PC1.x, PC2.x, PC3.x, PC4.x, PC5.x)
tmp <- tmp[complete.cases(tmp),] 
LM <- glm(IPF ~ HighPRS + telomere_gene + ShortTS + age.x + geneticSex.x + smoking + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x, family="binomial", data=data1)

out <- data.frame(summary(LM)$coefficients[-1,])
out <- out %>% mutate(variable = rownames(out))

out %>% write.table("category_IPF.susceptibility.assoc.tsv", sep="\t", quote=F, col.names = T, row.names = F)


