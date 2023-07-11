setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/")

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

baseline <- fread("~/gpfs1/01.IPFWGS/f2/01.basicQC/01.SampleQC/07.clinical/IPbaseline.txt")
ipf <- fread("../data/IPsample20220601.list")
ipf <- ipf %>% select(`CGM ID (二次)`, recruteday, sex, birthday)
colnames(ipf) <- c("ID", "DateOfRecruitment", "Sex", "DateOfBirth")
ipf <- ipf %>% mutate(DateOfRecruitment = as.Date(DateOfRecruitment, format="%d/%m/%Y"),
                      DateOfBirth = as.Date(DateOfBirth, format="%d/%m/%Y"),
                      Sex = ifelse(Sex == 0, "1|M", "2|F"))
ipf <- ipf %>% filter(!(ID %in% baseline$ID))

newbaseline <- bind_rows(baseline, ipf)

idlist <- openxlsx::read.xlsx("../data/IPFlist.xlsx")
idlist <- idlist %>% select(CGM_ID_2, KU_ID, death_date, cause_of_death)
idlist <- idlist %>% mutate(death_date = as.Date(as.numeric(death_date), origin = "1899-12-30"))
newbaseline <- newbaseline %>% left_join(idlist, by=c("ID"="CGM_ID_2"))
newbaseline <- newbaseline %>% mutate(KU_ID = case_when(ID == "PFKT0109" ~ 97804306,
                                                        ID == "PFKT0112" ~ 75853397,
                                                        ID == "PFKT0262" ~ 86853481,
                                                        ID == "PFKT0568" ~ 91504372,
                                                        ID == "PFKT0578" ~ 52423262,
                                                        ID == "PFKT1304" ~ 79663338,
                                                        ID == "PFKT1358" ~ 27383267,
                                                        ID == "PFKT2017" ~ 62884177,
                                                        ID == "PFKT2067" ~ 38044290,
                                                        ID == "PFKT2176" ~ 31533395,
                                                        ID == "PFKT2222" ~ 99163463,
                                                        ID == "PFKT2279" ~ 84434396,
                                                        ID == "PFKT2282" ~ 9053123,
                                                        ID == "PFKT2293" ~ 94724379,
                                                        ID == "PFKT2305" ~ 50794358,
                                                        ID == "PFKT2312_wd" ~ 48214350,
                                                        ID == "PFKT2334_wd" ~ 48844480,
                                                        ID == "PFKT2377_wd" ~ 93943247,
                                                        ID == "PFKT2405" ~ 53304364,
                                                        ID == "PFKT2438" ~ 4905189,
                                                        ID == "PFKT2442" ~ 31864306,
                                                        ID == "PFKT2513" ~ 61263506,
                                                        TRUE ~ KU_ID))
newbaseline <- newbaseline %>% mutate(DateOfDeath = case_when(is.na(as.Date(DateOfDeath)) ~ death_date,
                                                              TRUE ~ as.Date(DateOfDeath)))

tmp <- newbaseline %>% filter(!is.na(KU_ID) & (is.na(DateOfDiagnosis) | is.na(DateOfDeath)))　
tmp %>% write.xlsx("../data/NeedToExtract20220916.xlsx")

tmp <- read.xlsx("../data/NeedToExtract20220916_added.xlsx")
tmp <- tmp %>% mutate(DateOfRecruitment = as.Date(DateOfRecruitment, origin = "1899-12-30"),
                      DateOfBirth = as.Date(DateOfBirth, origin = "1899-12-30"),
                      DateOfDiagnosis = as.Date(as.numeric(DateOfDiagnosis), origin = "1899-12-30"),
                      DateOfDeath = as.Date(as.numeric(DateOfDeath), origin="1899-12-30"))
tmp1 <- newbaseline %>% filter(!(!is.na(KU_ID) & (is.na(DateOfDiagnosis) | is.na(DateOfDeath))))

newbaseline <- bind_rows(tmp, tmp1) %>% select(-death_date)

newbaseline %>% write.table("../data/IP_Diagnosis_Followup_20220916.tsv", sep="\t", quote=F, col.names = T, row.names = F)
