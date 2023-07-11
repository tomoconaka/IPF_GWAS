#Usage 
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)

before <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/GWAS/sIPF.SAIGEformat.txt.gz")
after <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/GWAS/sIPF.SAIGEformat.txt.gz")

SNPs <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/GWAS/sigloci", header=F)
before <- before %>% inner_join(SNPs, by=c("CHR"="V1", "POS"="V2"))
after  <- after  %>% inner_join(SNPs, by=c("CHR"="V1", "POS"="V2"))

tmp1 <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/GWAS/Kyoto_IPF.txt.gz")
tmp1 <- tmp1 %>% filter(Allele1 != "*" & Allele2 != "*")
dim(tmp1)#12117897

lambda_qc <- function(p){
  median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
}

lambda_qc(tmp1$p.value)#0.9457603



tmp2 <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/UKB_IPF.txt.gz")
tmp2 <- tmp2 %>% dplyr::select(p.value)

lambda_qc <- function(p){
  median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
}
lambda_qc(tmp1$p.value)
lambda_qc(tmp2$p.value)

png(file=paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/GWAS/GWAS_Kyoto.png"), width = 500, height = 500)
qq(tmp1$p.value)
dev.off()

filenames <- list.files(path="/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/",
                        pattern="*.SAIGEformat.txt.gz")

for(i in seq(1,6)){
  tmp <- fread(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/",filenames[i])) 
  lamda <- lambda_qc(tmp$p.value)
  title <- paste0(str_split(filenames[i], pattern="\\.", simplify = T)[,2], " ", str_split(filenames[i], pattern="\\.", simplify = T)[,4]) 
  title <- gsub("001", "", title)
  title <- gsub("01", "", title)
  title <- paste0(title, " λGC=",round(lamda, 2))
  png(file=paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/",filenames[i],".png"), width = 500, height = 500)
  qq(tmp$p.value, main = title)
  dev.off()
}

