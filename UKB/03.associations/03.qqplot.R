setwd("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/genebased")

lambda_qc <- function(p){
  median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
}

filenames <- list.files(pattern="*.SAIGEformat.txt.gz")

for(i in seq(1,6)){
  tmp <- fread(paste0(filenames[i])) 
  lamda <- lambda_qc(tmp$p.value)
  title <- paste0(str_split(filenames[i], pattern="\\.", simplify = T)[,2], " ", str_split(filenames[i], pattern="\\.", simplify = T)[,4]) 
  title <- gsub("001", "", title)
  title <- gsub("01", "", title)
  title <- paste0(title, " λGC=",round(lamda, 2))
  png(file=paste0(filenames[i],".png"), width = 500, height = 500)
  qq(tmp$p.value, main = title)
  dev.off()
}
