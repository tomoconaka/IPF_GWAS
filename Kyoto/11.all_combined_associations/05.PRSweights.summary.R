

meta <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/LeaveKyotoOut_meta_IPF.b38.txt.gz")
variants <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/IPF.extract")
variants <- variants %>% mutate(MarkerName = gsub("chr", "", paste0(V1,":",V2)))

meta <- meta %>% filter(MarkerName %in% variants$MarkerName)
meta <- meta %>% filter(!(MarkerName %in% c("5:169588475", "16:112241", "17:46253848")))

meta <- meta %>% select(CHR, POS, Allele1, Effect)
meta_JPN <- meta %>% mutate(EffectAllele = toupper(Allele1),
                            EffectWeight = Effect) %>% select(CHR, POS, EffectAllele, EffectWeight)

write.table(meta_JPN, file="/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Tables/PRSweight.JPN.tsv", 
            quote=F, col.names = T, row.names = F, sep="\t")

meta <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/LeaveUKBOut_meta_IPF.b38.txt.gz")
meta <- meta %>% filter(MarkerName %in% variants$MarkerName)
meta <- meta %>% filter(!(MarkerName %in% c("17:46258529")))
meta_EUR <- meta %>% mutate(EffectAllele = toupper(Allele1),
                            EffectWeight = Effect) %>% select(CHR, POS, EffectAllele, EffectWeight)

write.table(meta_EUR, file="/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Tables/PRSweight.EUR.tsv", 
            quote=F, col.names = T, row.names = F, sep="\t")


