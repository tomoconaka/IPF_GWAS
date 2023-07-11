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
##SBATCH --array=1-20,22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

module load ensembl-vep/104.3

SLURM_ARRAY_TASK_ID=X
awk -F"\t" '{print $4}' <(grep -v "#" $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.LoF.txt) | sort | uniq \
> $SCRATCH/annotation/LoF.gene.chr${SLURM_ARRAY_TASK_ID}.id

for gene in $(cat $SCRATCH/annotation/LoF.gene.chr${SLURM_ARRAY_TASK_ID}.id)
do
filter_vep -i $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.LoF.txt -o stdout -filter "Gene is ${gene}" | grep -v "#" | cut -f2 | sort | uniq | awk '($1 ~ /\-/){print $1}($1 !~ /\-/){a=substr($1, index($1, ":")+1); OFS=""; print $1,"-",a}' > $SCRATCH/annotation/tmp${SLURM_ARRAY_TASK_ID}
cat <(echo "${gene}") $SCRATCH/annotation/tmp${SLURM_ARRAY_TASK_ID} | paste -s | sed -e "s/\t/,/g" | sed -e "s/,/ /1" >> $SCRATCH/annotation/LoF.gene_list.chr${SLURM_ARRAY_TASK_ID}.txt
done

awk -F"\t" '{print $4}' <(grep -v "#" $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.CADD.txt) | sort | uniq \
> $SCRATCH/annotation/CADD.gene.chr${SLURM_ARRAY_TASK_ID}.id

for gene in $(cat $SCRATCH/annotation/CADD.gene.chr${SLURM_ARRAY_TASK_ID}.id)
do
filter_vep -i $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.CADD.txt -o stdout -filter "Gene is ${gene}" | grep -v "#" | cut -f2 | sort | uniq | awk '($1 ~ /\-/){print $1}($1 !~ /\-/){a=substr($1, index($1, ":")+1); OFS=""; print $1,"-",a}' > $SCRATCH/annotation/tmp${SLURM_ARRAY_TASK_ID}
cat <(echo "${gene}") $SCRATCH/annotation/tmp${SLURM_ARRAY_TASK_ID} | paste -s | sed -e "s/\t/,/g" | sed -e "s/,/ /1" >> $SCRATCH/annotation/CADD.gene_list.chr${SLURM_ARRAY_TASK_ID}.txt
rm $SCRATCH/annotation/tmp${SLURM_ARRAY_TASK_ID}
done
