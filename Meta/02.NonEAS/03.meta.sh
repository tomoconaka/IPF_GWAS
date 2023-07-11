#!/bin/bash


cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/
~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNP
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS GBMI_EAS.txt.gz
PROCESS Kyoto_IPF.txt.gz
PROCESS Allen_5way.txt.gz
PROCESS FinnGen_IPF.txt.gz
PROCESS UKB_IPF.txt.gz
PROCESS GBMI_AMR.txt.gz
PROCESS GBMI_AFR.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > ALL_meta_IPF.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header 

~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNP
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS GBMI_EAS.txt.gz
PROCESS Kyoto_ILD.txt.gz
PROCESS BBJv2_ILD.txt.gz
PROCESS Allen_5way.txt.gz
PROCESS FinnGen_ILD.txt.gz
PROCESS UKB_ILD.txt.gz
PROCESS GBMI_AMR.txt.gz
PROCESS GBMI_AFR.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > ALL_meta_ILD.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header

~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNP
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS Allen_IPF.b38_EAFfiexed.txt.gz
PROCESS FinnGen_IPF.txt.gz
PROCESS UKB_IPF.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > EUR_meta_IPF.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header


cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/
~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize

MARKER   MarkerName
WEIGHT   TotalSampleSize
ALLELE   Allele1 Allele2
EFFECT   Effect
FREQ     Freq1
STDERR   StdErr
PVAL     Pvalue

PROCESS EAS_meta_IPF.b38.txt.gz
PROCESS EUR_meta_IPF.b38.txt.gz 

LABEL TotalSampleSize as N

MARKER   SNP
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value
PROCESS GBMI_AMR.txt.gz
PROCESS GBMI_AFR.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > ALL_meta_IPF.percontinental.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header

cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/
~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNP
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS GBMI_EAS.txt.gz
PROCESS Allen_5way.txt.gz
PROCESS FinnGen_IPF.txt.gz
PROCESS UKB_IPF.txt.gz
PROCESS GBMI_AMR.txt.gz
PROCESS GBMI_AFR.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > LeaveKyotoOut_meta_IPF.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header

cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/
~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNP
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS GBMI_EAS.txt.gz
PROCESS Kyoto_IPF.txt.gz
PROCESS Allen_5way.txt.gz
PROCESS FinnGen_IPF.txt.gz
PROCESS GBMI_AMR.txt.gz
PROCESS GBMI_AFR.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > LeaveUKBOut_meta_IPF.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header


