
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
LoFadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased/IPF.LoF.0.001additive.b38.txt.gz")
LoFadd <- LoFadd %>% mutate(gene = str_split(MarkerName, pattern="\\.", simplify = T)[,1],
                            model = "LoF-additive")
LoFadd <- LoFadd %>% mutate(n_miss = str_count(Direction, "\\?"))
LoFadd <- LoFadd %>% filter(n_miss == 0)

LoFdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased/IPF.LoF.0.001dominant.b38.txt.gz")
LoFdom <- LoFdom %>% mutate(gene = str_split(MarkerName, pattern="\\.", simplify = T)[,1],
                            model = "LoF-dominant")
LoFdom <- LoFdom %>% mutate(n_miss = str_count(Direction, "\\?"))
LoFdom <- LoFdom %>% filter(n_miss == 0)

LoFrec <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased/IPF.LoF.0.01recessive.b38.txt.gz")
LoFrec <- LoFrec %>% mutate(gene = str_split(MarkerName, pattern="\\.", simplify = T)[,1],
                            model = "LoF-recessive")
LoFrec <- LoFrec %>% mutate(n_miss = str_count(Direction, "\\?"))
LoFrec <- LoFrec %>% filter(n_miss == 0)


CADDadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased/IPF.CADD20.0.001additive.b38.txt.gz")
CADDadd <- CADDadd %>% mutate(gene = str_split(MarkerName, pattern="\\.", simplify = T)[,1],
                            model = "LoF+CADD20-additive")
CADDadd <- CADDadd %>% mutate(n_miss = str_count(Direction, "\\?"))
CADDadd <- CADDadd %>% filter(n_miss == 0)

CADDdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased/IPF.CADD20.0.001dominant.b38.txt.gz")
CADDdom <- CADDdom %>% mutate(gene = str_split(MarkerName, pattern="\\.", simplify = T)[,1],
                            model = "LoF+CADD20-dominant")
CADDdom <- CADDdom %>% mutate(n_miss = str_count(Direction, "\\?"))
CADDdom <- CADDdom %>% filter(n_miss == 0)

CADDrec <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased/IPF.CADD20.0.01recessive.b38.txt.gz")
CADDrec <- CADDrec %>% mutate(gene = str_split(MarkerName, pattern="\\.", simplify = T)[,1],
                            model = "LoF+CADD20-recessive")
CADDrec <- CADDrec %>% mutate(n_miss = str_count(Direction, "\\?"))
CADDrec <- CADDrec %>% filter(n_miss == 0)

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
ALL <- ALL %>% select(gene, hgnc_symbol, chromosome_name, transcription_start_site, model, 
                      Effect, StdErr, Pvalue, Direction, HetPVal, TotalSampleSize)

ALL <- ALL %>% mutate(OR = exp(Effect),
               LL = exp(Effect + qnorm(0.025)*StdErr),
               UL = exp(Effect + qnorm(0.975)*StdErr))

ALL <- ALL %>% select(gene, hgnc_symbol, chromosome_name, transcription_start_site, model, 
                      OR, LL, UL, Pvalue, Direction, HetPVal, TotalSampleSize)


openxlsx::write.xlsx(ALL, file="/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Tables/Meta_genebased.xlsx")

ALL %>% filter(Pvalue < 0.05/dim(ALL)[1])

sig <- ALL %>% filter(Pvalue < 0.05/dim(ALL)[1])
sig <- sig %>% mutate(Gene = hgnc_symbol, 
                      Model = model,
                      OR_95CI = paste0(round(OR, 2), " [",round(LL, 2),"; ",round(UL, 2),"]"),
                      ensembleID = gene)

sig <- sig %>% select(Gene, ensembleID, chromosome_name, transcription_start_site, Model, OR_95CI, Pvalue,
                      HetPVal)
  
LoFadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/IPF.LoF.0.001additive.SAIGEformat.txt.gz")
LoFadd <- LoFadd %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-additive")
LoFadd <- LoFadd %>% filter(gene == "ENSG00000140694")
LoFdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/IPF.LoF.0.001dominant.SAIGEformat.txt.gz")
LoFdom <- LoFdom %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-dominant")
LoFdom <- LoFdom %>% filter(gene == "ENSG00000140694")

CADDadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/IPF.CADD20.0.001additive.SAIGEformat.txt.gz")
CADDadd <- CADDadd %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-additive")
CADDadd <- CADDadd %>% filter(gene %in% c("ENSG00000140694", "ENSG00000164362"))

CADDdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/IPF.CADD20.0.001dominant.SAIGEformat.txt.gz")
CADDdom <- CADDdom %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-dominant")
CADDdom <- CADDdom %>% filter(gene %in% c("ENSG00000140694", "ENSG00000164362"))

CADDrec <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/IPF.CADD20.0.01recessive.SAIGEformat.txt.gz")
CADDrec <- CADDrec %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-recessive")
CADDrec <- CADDrec %>% filter(gene == "ENSG00000176853")
Kyoto <- bind_rows(LoFadd, LoFdom, CADDadd, CADDdom, CADDrec)
Kyoto <- Kyoto %>% mutate(OR = exp(BETA),
                          LL = exp(BETA + qnorm(0.025)*SE),
                          UL = exp(BETA + qnorm(0.975)*SE))
Kyoto <- Kyoto %>% mutate(KyotoOR_95CI = paste0(round(OR, 2), " [",round(LL, 2),"; ",round(UL, 2),"]"),
                          Kyotopval =  p.value,
                          Kyoto_N = N, 
                          Kyoto_MAF = AF_Allele2)

Kyoto <- Kyoto %>% select(gene, model,KyotoOR_95CI, Kyotopval, Kyoto_MAF, Kyoto_N)

sig <- sig %>% inner_join(Kyoto, by=c("ensembleID"="gene",
                                      "Model"="model"))


LoFadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/genebased/IPF.LoF.0.001additive.SAIGEformat.txt.gz")
LoFadd <- LoFadd %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-additive")
LoFadd <- LoFadd %>% filter(gene == "ENSG00000140694")
LoFdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/genebased/IPF.LoF.0.001dominant.SAIGEformat.txt.gz")
LoFdom <- LoFdom %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                            model = "LoF-dominant")
LoFdom <- LoFdom %>% filter(gene == "ENSG00000140694")

CADDadd <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/genebased/IPF.CADD20.0.001additive.SAIGEformat.txt.gz")
CADDadd <- CADDadd %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-additive")
CADDadd <- CADDadd %>% filter(gene %in% c("ENSG00000140694", "ENSG00000164362"))

CADDdom <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/genebased/IPF.CADD20.0.001dominant.SAIGEformat.txt.gz")
CADDdom <- CADDdom %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-dominant")
CADDdom <- CADDdom %>% filter(gene %in% c("ENSG00000140694", "ENSG00000164362"))

CADDrec <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch/genebased/IPF.CADD20.0.01recessive.SAIGEformat.txt.gz")
CADDrec <- CADDrec %>% mutate(gene = str_split(SNPID, pattern="\\.", simplify = T)[,1],
                              model = "LoF+CADD20-recessive")
CADDrec <- CADDrec %>% filter(gene == "ENSG00000176853")
ukb <- bind_rows(LoFadd, LoFdom, CADDadd, CADDdom, CADDrec)
ukb <- ukb %>% mutate(OR = exp(BETA),
                      LL = exp(BETA + qnorm(0.025)*SE),
                      UL = exp(BETA + qnorm(0.975)*SE))
ukb <- ukb %>% mutate(UKBOR_95CI = paste0(round(OR, 2), " [",round(LL, 2),"; ",round(UL, 2),"]"),
                      UKBpval =  p.value,
                      UKB_N = N,
                      UKB_MAF = AF_Allele2)

ukb <- ukb %>% select(gene, model,UKBOR_95CI, UKBpval, UKB_MAF, UKB_N)
sig <- sig %>% inner_join(ukb, by=c("ensembleID"="gene",
                                      "Model"="model"))

openxlsx::write.xlsx(sig, file="../repo/IPF_GWAS/Tables/Table2.xlsx")


ALL <- readxl::read_excel("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Tables/Meta_genebased.xlsx")

candidates <- c("ABCA3", "DKC1", "KIF15", "MUC5B", "NAF1", "PARN", "RTEL1",
                "SFTPA1", "SFTPA2", "SFTPC", "TERC", "TERT",
                "GPR157", "MCL1", "MUC1", "THBS3", "RGS5", "PTPN14", "SERTAD2",
                "FAM13A", "SPDL1", "DSP", "HMGN4", "HLA-DRA","ZNF318","HSF2", "MAD1L1", "ZKSCAN1",
                "DEPTOR", "ZNF462", "STN1",  "DRAM1", "ATP11A", "IVD","GTPBP4",
                "AKAP13", "NPRL3", "SLC12A4", "ARL4D","KANSL1", "EVPL", "TYMS", "DPP9", "FUT3")


tmp <- ALL %>% filter(hgnc_symbol %in% candidates)

openxlsx::write.xlsx(tmp, file="/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Tables/SupplementaryTable5.xlsx")
length(unique(candidates))

tmp %>% filter(Pvalue < 0.05/44) %>% arrange(Pvalue)

