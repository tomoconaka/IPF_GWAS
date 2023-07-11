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
#SBATCH --array=2

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

module load ensembl-vep/104.3
mkdir -p $SCRATCH/annotation/tmp/

chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f1)
gene=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f2)
symbol=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SCRATCH/annotation/IPFgenes.list | cut -f3)
grep ${gene} $SCRATCH/annotation/${chr}.LoF.txt | awk '{print $1}' | sort | uniq > $SCRATCH/annotation/tmp/${gene}.LoF
zcat $SCRATCH/CADD.chr${chr}.txt.gz | grep "#" > $SCRATCH/annotation/tmp/${gene}tmp.txt
zcat $SCRATCH/CADD.chr${chr}.txt.gz | awk 'BEGIN{FS="\t"}($4 == "'${gene}'"){print $0}' >> $SCRATCH/annotation/tmp/${gene}tmp.txt
filter_vep -i $SCRATCH/annotation/tmp/${gene}tmp.txt -o $SCRATCH/annotation/tmp/${gene}Clinvar.txt -filter "CLIN_SIG is likely_pathogenic or CLIN_SIG is pathogenic" --force_overwrite
filter_vep -i $SCRATCH/annotation/tmp/${gene}tmp.txt -o $SCRATCH/annotation/tmp/${gene}common.txt -filter "MAX_AF > 0.01" --force_overwrite
grep -v "#" $SCRATCH/annotation/tmp/${gene}common.txt | awk -F"\t" '{print $1}' | sort | uniq > $SCRATCH/annotation/tmp/${gene}common
grep ${gene} $SCRATCH/annotation/${chr}.CADD.txt | awk '{print $1}' | sort | uniq > $SCRATCH/annotation/tmp/${gene}.CADD
grep ${gene} $SCRATCH/annotation/tmp/${gene}Clinvar.txt | awk '{print $1}' | sort | uniq > $SCRATCH/annotation/tmp/${gene}.Clinvar
comm -23 --nocheck-order $SCRATCH/annotation/tmp/${gene}.LoF $SCRATCH/annotation/tmp/${gene}common | awk '{OFS="\t"; m="LoF"; print m,$0}' > $SCRATCH/annotation/tmp/${symbol}.list
comm -23 --nocheck-order $SCRATCH/annotation/tmp/${gene}.CADD $SCRATCH/annotation/tmp/${gene}common | awk '{OFS="\t"; m="CADD"; print m,$0}' >> $SCRATCH/annotation/tmp/${symbol}.list
comm -23 --nocheck-order $SCRATCH/annotation/tmp/${gene}.Clinvar $SCRATCH/annotation/tmp/${gene}common | awk '{OFS="\t"; m="Clinvar"; print m,$0}' >> $SCRATCH/annotation/tmp/${symbol}.list
awk -F"\t" '{OFS="\t"; print $2}' $SCRATCH/annotation/tmp/${symbol}.list | sed -e "s/:/\t/g" | sed -e "s/-/\t/g" | sed -e "s/_/\t/g" | awk '{OFS="\t"; print $1,$2}'> $SCRATCH/annotation/tmp/${symbol}.range
#bcftools view ${SCRATCH}/FPF/${chr}.vcf.gz -R $SCRATCH/annotation/tmp/${symbol}.range -Oz -o ${SCRATCH}/FPF/${symbol}.vcf.gz
#bcftools view ${SCRATCH}/IPF/${chr}.vcf.gz -R $SCRATCH/annotation/tmp/${symbol}.range -Oz -o ${SCRATCH}/IPF/${symbol}.vcf.gz
