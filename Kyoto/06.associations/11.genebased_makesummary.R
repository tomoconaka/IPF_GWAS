
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
LoFadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/genebased/IPF.LoF.0.001additive.SAIGEformat.txt.gz")
LoFadd <- LoFadd %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-additive")

LoFdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/genebased/IPF.LoF.0.001dominant.SAIGEformat.txt.gz")
LoFdom <- LoFdom %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-dominant")

LoFrec <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/genebased/IPF.LoF.0.01recessive.SAIGEformat.txt.gz")
LoFrec <- LoFrec %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-recessive")


CADDadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/genebased/IPF.CADD20.0.001additive.SAIGEformat.txt.gz")
CADDadd <- CADDadd %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-additive")


CADDdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/genebased/IPF.CADD20.0.001dominant.SAIGEformat.txt.gz")
CADDdom <- CADDdom %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-dominant")

CADDrec <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/genebased/IPF.CADD20.0.01recessive.SAIGEformat.txt.gz")
CADDrec <- CADDrec %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-recessive")


lambda_qc <- function(p){
  median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
}

png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/ExtendedDataFig.KyotoQQLoFadd.png")), width=500, height=500)
qqman::qq(LoFadd$p.value,
          main = paste0("Kyoto ",unique(LoFadd$model)," model (λGC = ",round(lambda_qc(LoFadd$p.value), 2),")"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)

dev.off()

png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/ExtendedDataFig.KyotoQQLoFdom.png")), width=500, height=500)
qqman::qq(LoFdom$p.value,
          main = paste0("Kyoto ",unique(LoFdom$model)," model (λGC = ",round(lambda_qc(LoFdom$p.value), 2),")"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)

dev.off()

png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/ExtendedDataFig.KyotoQQLoFrec.png")), width=500, height=500)
qqman::qq(LoFrec$p.value,
          main = paste0("Kyoto ",unique(LoFrec$model)," model (λGC = ",round(lambda_qc(LoFrec$p.value), 2),")"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)

dev.off()

#############
png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/ExtendedDataFig.KyotoQQCADD20add.png")), width=500, height=500)
qqman::qq(CADDadd$p.value,
          main = paste0("Kyoto ",unique(CADDadd$model)," model (λGC = ",round(lambda_qc(CADDadd$p.value), 2),")"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)

dev.off()

png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/ExtendedDataFig.KyotoQQCADD20dom.png")), width=500, height=500)
qqman::qq(CADDdom$p.value,
          main = paste0("Kyoto ",unique(CADDdom$model)," model (λGC = ",round(lambda_qc(CADDdom$p.value), 2),")"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)

dev.off()

png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/ExtendedDataFig.KyotoQQCADD20rec.png")), width=500, height=500)
qqman::qq(CADDrec$p.value,
          main = paste0("Kyoto ",unique(CADDrec$model)," model (λGC = ",round(lambda_qc(CADDrec$p.value), 2),")"), xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)

dev.off()


ALL <- bind_rows(LoFadd, LoFdom, LoFrec, CADDadd, CADDdom, CADDrec)


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(ALL$gene))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("chromosome_name","ensembl_gene_id", "ensembl_transcript_id",'hgnc_symbol', "transcription_start_site", "transcript_mane_select"),
                filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)
G_list1 <- G_list %>% arrange(desc(transcript_mane_select), transcription_start_site)
G_list1 <- G_list1[!duplicated(G_list1$ensembl_gene_id),]

ALL <- ALL %>% left_join(G_list1, by=c("gene"="ensembl_gene_id"))
ALL <- ALL %>% dplyr::select(gene, hgnc_symbol, chromosome_name, transcription_start_site, AF_Allele2, BETA, SE, p.value, model)

ALL <- ALL %>% rename(EnsembleID = gene, 
                      CHR = chromosome_name,
                      TSS = transcription_start_site,
                      AlleleFrequency = AF_Allele2)
ALL <- ALL %>% mutate(OR = exp(BETA), 
                      LL = exp(BETA + qnorm(0.025)*SE),
                      UL = exp(BETA + qnorm(0.975)*SE))

ALL <- ALL %>% mutate(EnsembleID, CHR, TSS, AlleleFrequency, OR, LL, UL, p.value, model)

openxlsx::write.xlsx(ALL, file="../repo/IPF_GWAS/Tables/Kyoto_genebased.xlsx")
