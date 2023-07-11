#!/bin/bash


~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNP
ALLELE   Allele2 Allele1
WEIGHTLABEL     N
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/GBMI_EAS.txt.gz
PROCESS /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/Kyoto_IPF.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/EAS_meta_IPF.b38.txt.gz

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

PROCESS /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/BBJv2_ILD.txt.gz
PROCESS /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/GBMI_EAS.txt.gz
PROCESS /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/Kyoto_ILD.txt.gz

ANALYZE RANDOM
###

paste <(echo "CHR") <(echo "POS") <(head -1 METAANALYSIS1.TBL) -d "\t" > meta.header
cat meta.header > meta.result
paste <(awk '{print $1}' METAANALYSIS1.TBL | tail -n+2 | sed -e "s/\:/\t/g") <(awk '{OFS="\t"; print $0}' METAANALYSIS1.TBL | tail -n+2) -d "\t" | sed -e "s/X/23/g" | awk '$1 ~ /^[0-9]+$/' |  sort -k1,1n -k2,2n >> meta.result
cat meta.result | bgzip -c > /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/GBMI_Kyoto_BBJv2_ILD.b38.txt.gz

rm METAANALYSIS1.TBL*
rm meta.header 

