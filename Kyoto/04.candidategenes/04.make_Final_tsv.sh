#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=vep
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-13

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

#module load ensembl-vep/104.3
mkdir -p $SCRATCH/annotation/tmp/


for pheno in FPF IPF AGP3000
do
for SLURM_ARRAY_TASK_ID in 1
do
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f1)
gene=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f2)
symbol=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f3)
rm ${SCRATCH}/${pheno}/${pheno}.Final.tsv
zcat ${SCRATCH}/${pheno}/${symbol}.vcf.gz | grep -v "##" | head -1 > ${SCRATCH}/${pheno}/${symbol}.final.tsv
zcat ${SCRATCH}/${pheno}/${symbol}.vcf.gz | grep -v "#" | grep -e "\/1" -e "0|1" -e "1|1" >> ${SCRATCH}/${pheno}/${symbol}.final.tsv
cat ${SCRATCH}/${pheno}/${symbol}.final.tsv > ${SCRATCH}/${pheno}/${pheno}.Final.tsv
done
done

for pheno in FPF IPF AGP3000
do
for SLURM_ARRAY_TASK_ID in {2..13}
do
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f1)
gene=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f2)
symbol=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f3)
zcat ${SCRATCH}/${pheno}/${symbol}.vcf.gz | grep -v "##" | head -1 > ${SCRATCH}/${pheno}/${symbol}.final.tsv
zcat ${SCRATCH}/${pheno}/${symbol}.vcf.gz | grep -v "#" | grep -e "\/1" -e "0|1" -e "1|1" >> ${SCRATCH}/${pheno}/${symbol}.final.tsv
tail -n+2 ${SCRATCH}/${pheno}/${symbol}.final.tsv >> ${SCRATCH}/${pheno}/${pheno}.Final.tsv
done
done
