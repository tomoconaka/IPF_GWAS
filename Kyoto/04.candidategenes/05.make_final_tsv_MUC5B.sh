#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=vep
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-13

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

#module load ensembl-vep/104.3
mkdir -p $SCRATCH/annotation/tmp/

for pheno in FPF IPF AGP3000
do
bcftools view ${SCRATCH}/${pheno}/11.vcf.gz -r chr11:1219991 -Oz -o ${SCRATCH}/${pheno}/MUC5B.vcf.gz
zcat ${SCRATCH}/${pheno}/MUC5B.vcf.gz | grep -v "##" > ${SCRATCH}/${pheno}/MUC5B.final.tsv
done
#awk -F"\t" 'FNR==NR {m[$1":"$2]=$0; next} ($2 in m) {OFS="\t"; print "'${symbol}'",$7,$14,m[$2]}' \
#${SCRATCH}/FPF/${symbol}.final.tsv $SCRATCH/annotation/${chr}.LoF.txt $SCRATCH/annotation/${chr}.CADD.txt | sort -r | uniq -f3 >> ${SCRATCH}/FPF/Final.tsv

#zcat ${SCRATCH}/IPF/${symbol}.vcf.gz | grep -v "##" | head -1 > ${SCRATCH}/IPF/${symbol}.final.tsv
#zcat ${SCRATCH}/IPF/${symbol}.vcf.gz | grep -v "#" | grep -e "\/1" -e "0|1" -e "1|1" >> ${SCRATCH}/IPF/${symbol}.final.tsv
#awk -F"\t" 'FNR==NR {m[$1":"$2]=$0; next} ($2 in m) {OFS="\t"; print "'${symbol}'",$7,$14,m[$2]}' \
#${SCRATCH}/IPF/${symbol}.final.tsv $SCRATCH/annotation/${chr}.LoF.txt $SCRATCH/annotation/${chr}.CADD.txt | sort -r | uniq -f3 >> ${SCRATCH}/IPF/Final.tsv

#zcat ${SCRATCH}/AGP3000/${symbol}.vcf.gz | grep -v "##" | head -1 > ${SCRATCH}/AGP3000/${symbol}.final.tsv
#zcat ${SCRATCH}/AGP3000/${symbol}.vcf.gz | grep -v "#" | grep -e "\/1" -e "0|1" -e "1|1" >> ${SCRATCH}/AGP3000/${symbol}.final.tsv
#awk -F"\t" 'FNR==NR {m[$1":"$2]=$0; next} ($2 in m) {OFS="\t"; print "'${symbol}'",$7,$14,m[$2]}' \
#${SCRATCH}/IPF/${symbol}.final.tsv $SCRATCH/annotation/${chr}.LoF.txt $SCRATCH/annotation/${chr}.CADD.txt | sort -r | uniq -f3 >> ${SCRATCH}/AGP3000/Final.tsv

