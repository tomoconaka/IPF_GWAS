#!/bin/bash
#SBATCH --time=150:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=candidategenes
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

SLURM_ARRAY_TASK_ID=X
mkdir -p ${SCRATCH}/FPF
bcftools view -S ${SCRATCH}/../data/FPFsample --force-samples ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.2.vcf.gz -Oz -o ${SCRATCH}/FPF/${SLURM_ARRAY_TASK_ID}.vcf.gz
tabix -p vcf ${SCRATCH}/FPF/${SLURM_ARRAY_TASK_ID}.vcf.gz
mkdir -p ${SCRATCH}/IPF
bcftools view -S ${SCRATCH}/../data/IPFsample --force-samples ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.2.vcf.gz -Oz -o ${SCRATCH}/IPF/${SLURM_ARRAY_TASK_ID}.vcf.gz
tabix -p vcf ${SCRATCH}/IPF/${SLURM_ARRAY_TASK_ID}.vcf.gz
mkdir -p ${SCRATCH}/AGP3000
bcftools view -S ${SCRATCH}/../data/AGP3000sample --force-samples ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.2.vcf.gz -Oz -o ${SCRATCH}/AGP3000/${SLURM_ARRAY_TASK_ID}.vcf.gz
tabix -p vcf ${SCRATCH}/AGP3000/${SLURM_ARRAY_TASK_ID}.vcf.gz

