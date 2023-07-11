#!/bin/bash

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch

cd ${SCRATCH}
awk '(NR > 1){print $2}' ukb23149_450k_OQFE.variant_ID_mappings.txt | awk -F':' '($1 != 23 && $1 != 24){OFS="\t"; print "chr"$1, $2, $2+length($3)-1, $3"/"$4, "+", "chr"$1"_"$2"_"$3"_"$4} ($1 == 23) {OFS="\t"; print "chrX", $2, $2+length($3)-1, $3"/"$4, "+", "chrX_"$2"_"$3"_"$4} ($1 == 24){OFS="\t"; print "chrY", $2, $2+length($3)-1, $3"/"$4, "+", "chrY_"$2"_"$3"_"$4}' | bgzip -c > ukb23149_450k_OQFE.variant_ID.input.txt.gz

for chr in {1..22}
do
awk '(NR > 1){print $2}' ukb23149_450k_OQFE.variant_ID_mappings.txt | \
awk -F':' '($1 == '"$chr"'){OFS="\t"; print "chr"$1, $2, $2+length($3)-1, $3"/"$4, "+", "chr"$1"_"$2"_"$3"_"$4}' | bgzip -c > chr${chr}.ukb23149_450k_OQFE.variant_ID.input.txt.gz
done 

awk '(NR > 1){print $2}' ukb23149_450k_OQFE.variant_ID_mappings.txt | awk -F':' ' ($1 == 23) {OFS="\t"; print "chrX", $2, $2+length($3)-1, $3"/"$4, "+", "chrX_"$2"_"$3"_"$4} ($1 == 24){OFS="\t"; print "chrY", $2, $2+length($3)-1, $3"/"$4, "+", "chrY_"$2"_"$3"_"$4}' | bgzip -c > chrX.ukb23149_450k_OQFE.variant_ID.input.txt.gz

tail -n+2 ${SCRATCH}/ukb23149_450k_OQFE.variant_ID_mappings.txt | cut -f2 | awk -F':' '{OFS="\t"; print $0, $1, $2, $3, $4}' | \
awk '{OFS="\t"; print "chr"$2, $3, "chr"$2"_"$3"_"$4"_"$5,$4, $5, 100, ".", ".", "GT:AD:DP:GQ:PL", "0/0:1,0:1:3:0,3,15"}' \
>> $SCRATCH/ukb23149_450k_OQFE.variant.vcf
bgzip $SCRATCH/ukb23149_450k_OQFE.variant.vcf
tabix -p vcf $SCRATCH/ukb23149_450k_OQFE.variant.vcf.gz
