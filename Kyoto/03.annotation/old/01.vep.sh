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
##SBATCH --array=1

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

module load ensembl-vep/104.3

SLURM_ARRAY_TASK_ID=X
vep -i ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.2.vcf.gz \
--plugin dbNSFP,/home/tnakanishi/gpfs1/data/vep/dbNSFP4/dbNSFP4.3a_grch38.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
-everything \
--force_overwrite \
--buffer_size 10000 \
--offline \
--fork 15 \
--assembly GRCh38 \
--dir_plugins /home/tnakanishi/gpfs1/data/vep/VEP_plugins \
--dir_cache /home/tnakanishi/gpfs1/data/vep/cache/ \
--cache -o $SCRATCH/finalAnnot.chr${SLURM_ARRAY_TASK_ID}.txt
