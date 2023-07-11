#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=telseq
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p ALL
#SBATCH --array=3227-3228,3244-3255,3357-3366

cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/Telseq
source /home/tnakanishi/anaconda3/bin/activate telseq
ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list)
CRAM=$(grep ${ID} /home/nagasaki/workspace/pipeline/dist/v3/cramlist.v3.abs.txt)
samtools view -b ${CRAM} -T /mnt/gpfs1/home/nagasaki/workspace/pipeline/data/hs38DH.fa -o bams/${ID}.bam
telseq -r 150 bams/${ID}.bam > ${ID}.result

