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
##SBATCH --array=1-22

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch

cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Kyoto/06.associations

./regenie \
--step 1 \
--bed ${SCRATCH}/grm_final \
--phenoFile ${SCRATCH}/pheno.txt \
--phenoColList IPF,ILD,FPF,sIPF \
--covarFile ${SCRATCH}/pheno.txt \
--covarColList age,PC1,PC2,PC3,PC4,PC5 \
--catCovarList geneticSex,platform \
--bt \
--bsize 1000 \
--threads 20 \
--lowmem \
--lowmem-prefix ${SCRATCH}/GWAS/tmpdir/pheno_tmp_preds \
--out ${SCRATCH}/GWAS/pheno_step1 

