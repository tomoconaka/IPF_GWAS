#!/bin/bash
#SBATCH --time=100:00:00
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
##SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch

module load ensembl-vep/104.3

#SLURM_ARRAY_TASK_ID=X
filter_vep -i $SCRATCH/CADD.ukb23149_450k_OQFE.variant_ID.txt.gz \
-filter "CADD_PHRED > 20" \
-o $SCRATCH/annotation/ukb23149_450k_OQFE.variant_ID.CADD.txt --force_overwrite

