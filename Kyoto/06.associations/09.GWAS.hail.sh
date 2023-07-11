#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=hail
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH -p ALL
##SBATCH --array=1-22

source ~/anaconda3/bin/activate hail
export PYSPARK_SUBMIT_ARGS="--driver-memory 400G pyspark-shell"
export JAVA_OPTS="-Xss128k -Xms256m -Xmx512m -XX:PermSize=64m -XX:MaxPermSize=1024m"

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/GWAS/

python GWASplot.py --input ${SCRATCH}/${1} --outdir /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/  --pheno ${1} --qqplot
