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

cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/ILD_GWAS/06.associations

./regenie \
--step 2 \
--bed ${SCRATCH}/chr${SLURM_ARRAY_TASK_ID} \
--phenoFile ${SCRATCH}/pheno.txt \
--phenoColList IPF,ILD,FPF,sIPF \
--covarFile ${SCRATCH}/pheno.txt \
--covarColList age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--catCovarList geneticSex \
--bt \
--bsize 400 \
--threads 20 \
--firth --approx --firth-se \
--pThresh 0.01 \
--pred ${SCRATCH}/GWAS/pheno_step1_pred.list \
--out ${SCRATCH}/GWAS/chr${SLURM_ARRAY_TASK_ID}

