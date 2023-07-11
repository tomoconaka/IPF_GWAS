#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sampleQC
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

vcftools --gzvcf $SCRATCH/allchr.sampleQC.vcf.gz --missing-indv
vcftools --gzvcf $SCRATCH/allchr.sampleQC.vcf.gz --depth

#cat $SCRATCH/step1.tokeep.IP.rev.sample | grep -v "PFKT0797" | grep -v "PFKT0893" > $SCRATCH/step1.tokeep_rev.IP.sample

#cat $SCRATCH/step1.tokeep_rev.IP.sample $SCRATCH/step1.tokeep.AGP3000.sample > $SCRATCH/SampleQC.tokeep.All.sample

