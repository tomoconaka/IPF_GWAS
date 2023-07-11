#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=pca
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --array=1-22

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch


for i in {1..22}
do
echo ${SCRATCH}/grm_pass1/${i}_LDp >> ${SCRATCH}/grm_pass1/merge.list
done

plink \
--merge-list ${SCRATCH}/grm_pass1/merge.list \
--biallelic-only \
--list-duplicate-vars suppress-first \
--snps-only just-acgt \
--memory  240000 \
--threads 4 \
--make-bed \
--out ${SCRATCH}/grm_final

awk '{print $2,$2,$3,$4,$5,$6}' ${SCRATCH}/grm_final.fam > ${SCRATCH}/tmp
mv ${SCRATCH}/tmp ${SCRATCH}/grm_final.fam

#cat ${SCRATCH}/grm_final.fam | sed -e "s/C AC/C_AC/g" > ${SCRATCH}/tmp
#awk '($2 ~ /^C/){print $1,$0}($2 !~ /^C/){print $0}' ${SCRATCH}/tmp > ${SCRATCH}/grm_final.fam

#cat ${SCRATCH}/grm_final.fam | sed -e "s/ wd/_wd/g" > ${SCRATCH}/tmp
#awk '($2 ~ /wd/){print $1,$0}($2 !~ /wd/){print $0}' ${SCRATCH}/tmp > ${SCRATCH}/grm_final.fam

