#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sampleQC
##SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch
vcftools --gzvcf ${DATADIR}/all.VQSR3.chr${SLURM_ARRAY_TASK_ID}.vcf.gz --keep <(cat $SCRATCH/step1.tokeep.AGP3000.sample $SCRATCH/step1.tokeep.IP.rev.sample) --remove-indels --remove-filtered-all --maf 0.01 --max-missing 0.99 --recode --stdout | bgzip -c > $SCRATCH/chr${SLURM_ARRAY_TASK_ID}.sampleQC.vcf.gz
tabix -p vcf $SCRATCH/chr${SLURM_ARRAY_TASK_ID}.sampleQC.vcf.gz

#vcftools --gzvcf $SCRATCH/chr1-22.sampleQC.vcf.gz --missing-indv

#vcftools --gzvcf $SCRATCH/chr1-22.sampleQC.vcf.gz --depth

