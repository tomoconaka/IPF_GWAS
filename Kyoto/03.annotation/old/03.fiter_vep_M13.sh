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
#SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

module load ensembl-vep/104.3

filter_vep -i $SCRATCH/finalAnnot.chr${SLURM_ARRAY_TASK_ID}.txt \
-filter "Consequence is stop_gained or Consequence is start_lost or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is stop_lost or Consequence is frameshift_variant" \
-o $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.M1.txt --force_overwrite

filter_vep -i $SCRATCH/finalAnnot.chr${SLURM_ARRAY_TASK_ID}.txt \
-filter "Consequence is missense_variant and SIFT_pred is D and Polyphen2_HDIV_pred is D and Polyphen2_HVAR_pred is D and LRT_pred is D and (MutationTaster_pred is D or MutationTaster_pred is A)" \
-o $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.M3.txt --force_overwrite
