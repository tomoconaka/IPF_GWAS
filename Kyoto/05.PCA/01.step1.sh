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

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

mkdir -p ${SCRATCH}/grm_pass1/
plink2 --vcf ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.final.vcf.gz \
--maf 0.05 \
--geno 0.05 \
--snps-only just-acgt \
--indep-pairwise 50 5 0.05 \
--exclude range /home/tnakanishi/gpfs1/data/LdRegion-AbecasisHg38.txt \
--max-alleles 2 \
--memory 24000 \
--threads 8 \
--out ${SCRATCH}/grm_pass1/${SLURM_ARRAY_TASK_ID}_LDp

plink2 \
--vcf ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.final.vcf.gz \
--extract ${SCRATCH}/grm_pass1/${SLURM_ARRAY_TASK_ID}_LDp.prune.in \
--keep-allele-order \
--make-bed \
--out ${SCRATCH}/grm_pass1/${SLURM_ARRAY_TASK_ID}_LDp

