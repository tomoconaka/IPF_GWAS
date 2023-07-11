#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=vep
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p NGS8,NGS9,DGX1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=19-20,22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

module load ensembl-vep/104.3

bcftools view -R ~/gpfs1/data/1000G/LCR/b38.bed \
$SCRATCH/ukb23149_450k_OQFE.variant.vcf.gz \
-Oz -o $SCRATCH/ukb23149_450k_OQFE.variant.LCR.vcf.gz

for MAF in 0.001 0.01
do
filter_vep -i $SCRATCH/CADD.ukb23149_450k_OQFE.variant_ID.txt.gz \ 
-o $SCRATCH/${MAF}.txt -filter "MAX_AF > ${MAF}" --force_overwrite
done
