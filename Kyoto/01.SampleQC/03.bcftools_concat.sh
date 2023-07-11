#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sampleQC
##SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

bcftools concat -f $SCRATCH/chr1-22.sampleQC.vcf.list \
-Oz -o $SCRATCH/allchr.sampleQC.vcf.gz
