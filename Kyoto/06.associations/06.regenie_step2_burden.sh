#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=regenie_BBJ
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH -p ALL
#SBATCH --array=1-22

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch

cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Kyoto/06.associations

for SLURM_ARRAY_TASK_ID in {1..22} X PAR
do
for TEST in dominant additive
do
MAF=0.001

#SLURM_ARRAY_TASK_ID=PAR
./regenie \
--step 2 \
--minMAC 1 \
--test ${TEST} \
--bed ${SCRATCH}/${SLURM_ARRAY_TASK_ID} \
--aaf-bins ${MAF} \
--build-mask 'comphet' \
--phenoFile ${SCRATCH}/pheno.txt \
--phenoColList IPF,ILD,FPF,sIPF \
--covarFile ${SCRATCH}/pheno.txt \
--covarColList age,PC1,PC2,PC3,PC4,PC5 \
--catCovarList geneticSex,platform \
--mask-def ../../../../data/maskfile \
--exclude-sets ${SCRATCH}/genebased/chr${SLURM_ARRAY_TASK_ID}_exclude_sets_${MAF} \
--set-list ${SCRATCH}/genebased/set_list.chr${SLURM_ARRAY_TASK_ID}.txt \
--anno-file $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.txt \
--bt \
--bsize 400 \
--threads 20 \
--firth --approx --firth-se \
--pThresh 0.01 \
--pred ${SCRATCH}/GWAS/pheno_step1_pred.list \
--out ${SCRATCH}/genebased/chr${SLURM_ARRAY_TASK_ID}_${TEST}

done

for TEST in recessive
do
MAF=0.01

#SLURM_ARRAY_TASK_ID=PAR
./regenie \
--step 2 \
--minMAC 1 \
--test ${TEST} \
--bed ${SCRATCH}/${SLURM_ARRAY_TASK_ID} \
--aaf-bins ${MAF} \
--build-mask 'comphet' \
--phenoFile ${SCRATCH}/pheno.txt \
--phenoColList IPF,ILD,FPF,sIPF \
--covarFile ${SCRATCH}/pheno.txt \
--covarColList age,PC1,PC2,PC3,PC4,PC5 \
--catCovarList geneticSex,platform \
--mask-def ../../../../data/maskfile \
--exclude-sets ${SCRATCH}/genebased/chr${SLURM_ARRAY_TASK_ID}_exclude_sets_${MAF} \
--set-list ${SCRATCH}/genebased/set_list.chr${SLURM_ARRAY_TASK_ID}.txt \
--anno-file $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.txt \
--bt \
--bsize 400 \
--threads 20 \
--firth --approx --firth-se \
--pThresh 0.01 \
--pred ${SCRATCH}/GWAS/pheno_step1_pred.list \
--out ${SCRATCH}/genebased/chr${SLURM_ARRAY_TASK_ID}_${TEST}

done
