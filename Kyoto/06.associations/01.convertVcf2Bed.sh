#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH -p ALL
#SBATCH --array=1-22

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch

#SLURM_ARRAY_TASK_ID=PAR

plink --make-bed \
  --out ${SCRATCH}/${SLURM_ARRAY_TASK_ID} \
  --vcf ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.final.vcf.bgz

cat ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.fam | sed -e "s/C AC/C_AC/g" > ${SCRATCH}/tmp
awk '($1 ~ /^C/){print $1,$0}($1 !~ /^C/){print $0}' ${SCRATCH}/tmp > ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.fam

cat ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.fam | sed -e "s/ wd/_wd/g" > ${SCRATCH}/tmp
awk '($1 ~ /wd/){print $1,$0}($1 !~ /wd/){print $0}' ${SCRATCH}/tmp > ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.fam

plink --bfile ${SCRATCH}/${SLURM_ARRAY_TASK_ID} --maf 0.001 --geno 0.03 --make-bed --out ${SCRATCH}/${SLURM_ARRAY_TASK_ID}_for_lamba
