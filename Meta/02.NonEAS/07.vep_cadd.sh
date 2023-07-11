#!/bin/bash
#SBATCH --time=150:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=30
#SBATCH --job-name=vep
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p NGS8,NGS9,DGX1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-22

#DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch

#SLURM_ARRAY_TASK_ID=PAR
module load ensembl-vep/104.3

vep -i ${SCRATCH}/significant_loci.txt.gz \
--plugin CADD,/home/tnakanishi/gpfs1/data/vep/CADD/whole_genome_SNVs.tsv.gz \
--force_overwrite \
--buffer_size 10000 \
--offline \
--fork 20 \
--assembly GRCh38 \
-everything \
--dir_plugins /home/tnakanishi/gpfs1/data/vep/VEP_plugins \
--dir_cache /home/tnakanishi/gpfs1/data/vep/cache/ \
--tab --compress_output bgzip \
--cache -o $SCRATCH/significant_loci.annotated.txt.gz
