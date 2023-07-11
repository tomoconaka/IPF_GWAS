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

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

plink --bfile ${SCRATCH}/chr${SLURM_ARRAY_TASK_ID} \
--extract 'range' <(awk '{print $0, "L"NR}' ~/gpfs1/data/1000G/LCR/b38.bed) \
--write-snplist \
--out ${SCRATCH}/chr${SLURM_ARRAY_TASK_ID}.toexclude
