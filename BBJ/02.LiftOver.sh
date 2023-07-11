#!/bin/bash

SCRATCH=../../../BBJ_scratch/

zcat ${SCRATCH}/ILD.SAIGEformat.txt.gz | tail -n+2 | awk '{OFS="\t" ; print "chr"$1,int($2),int($2+1),$3}' | sed -e "s/chr23/chrX/g" > ${SCRATCH}/BBJ_ILD.b37.bed
liftOver ${SCRATCH}/BBJ_ILD.b37.bed ~/gpfs1/data/LiftOverchains/hg19ToHg38.over.chain ${SCRATCH}/BBJ_ILD.b38.bed ${SCRATCH}/unMapped

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 N" > ../../../Meta_scratch/BBJv2_ILD.txt
awk -F " " 'FNR==NR { m[$4]=$1; n[$4] = $2; next } \
($3 in m && $6 >= 0.001 && $6 <= 0.999 && $7 >= 0.3)  { OFS=" "; print m[$3]":"n[$3],$4,$5,$10,$11,$16,$6, $8}' ../../../Meta_scratch/BBJ_ILD.b38.bed | <(zcat $SCRATCH/ILD.SAIGEformat.txt.gz) |\
sed -e "s/chr//g" | sed -e "s/X:/23:/g" >> ../../../Meta_scratch/BBJv2_ILD.txt

