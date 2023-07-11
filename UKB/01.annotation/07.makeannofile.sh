#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=vep
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p NGS8,NGS9
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=1-22

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/UKB_scratch

module load ensembl-vep/104.3

SLURM_ARRAY_TASK_ID=1
cat $SCRATCH/annotation/${SLURM_ARRAY_TASK_ID}.LoF.txt | grep -v "#" | cut -f1,4 | grep -v "-" | sort | uniq > $SCRATCH/annotation/LoFtmp${SLURM_ARRAY_TASK_ID}
awk '{OFS="\t"; print $0,"LoF"}' $SCRATCH/annotation/LoFtmp${SLURM_ARRAY_TASK_ID} > $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.txt
cat $SCRATCH/annotation/ukb23149_450k_OQFE.variant_ID.CADDmissense.txt | grep -v "#" | awk '$2 ~ /chr'"${SLURM_ARRAY_TASK_ID}"':/' | cut -f1,4 | grep -v "-" | sort | uniq > $SCRATCH/annotation/CADDmissensetmp${SLURM_ARRAY_TASK_ID}
comm -23 --nocheck-order $SCRATCH/annotation/CADDmissensetmp${SLURM_ARRAY_TASK_ID} $SCRATCH/annotation/LoFtmp${SLURM_ARRAY_TASK_ID} | awk '{OFS="\t"; print $0,"CADDmissense"}' >> $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.txt
cat $SCRATCH/annotation/ukb23149_450k_OQFE.variant_ID.CADD.txt | grep -v "#" | awk '$2 ~ /chr'"${SLURM_ARRAY_TASK_ID}"':/' | cut -f1,4 | grep -v "-" | sort | uniq > $SCRATCH/annotation/CADDtmp${SLURM_ARRAY_TASK_ID}
comm -23 --nocheck-order $SCRATCH/annotation/CADDtmp${SLURM_ARRAY_TASK_ID} $SCRATCH/annotation/CADDmissensetmp${SLURM_ARRAY_TASK_ID} | awk '{OFS="\t"; print $0,"CADD"}' >> $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.txt

awk '{a[$1,$2]++}!(a[$1,$2]-1)' $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.txt | sed -e "s/\\_/\\:/g" > $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.rev.txt
cut -f2 $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.rev.txt | sort | uniq > ${SCRATCH}/genebased/settmp${SLURM_ARRAY_TASK_ID}
for i in $(cat ${SCRATCH}/genebased/settmp${SLURM_ARRAY_TASK_ID})
do
awk '$2 == "'$i'"{print $1}' $SCRATCH/annotation/annofile.chr${SLURM_ARRAY_TASK_ID}.rev.txt | tr "\n" "," | sed -e "s/,\$//" > ${SCRATCH}/genebased/set1tmp${SLURM_ARRAY_TASK_ID}
awk -F ":" '{OFS=" "; print "'$i'", $1, $2, $0}' ${SCRATCH}/genebased/set1tmp${SLURM_ARRAY_TASK_ID} | sed -e "s/chr//1" >> ${SCRATCH}/genebased/set_list.chr${SLURM_ARRAY_TASK_ID}.txt
done

for MAF in 0.001 0.01
do
grep -v "#" $SCRATCH/${MAF}.txt | awk '$2 ~ /chr'"${SLURM_ARRAY_TASK_ID}"':/' | cut -f1 | sort | uniq > $SCRATCH/${SLURM_ARRAY_TASK_ID}.tmp1
zcat $SCRATCH/ukb23149_450k_OQFE.variant.LCR.vcf.gz | grep -v "#" | awk '$1 == "chr'"${SLURM_ARRAY_TASK_ID}"'"' | cut -f3 > $SCRATCH/${SLURM_ARRAY_TASK_ID}.tmp
cat $SCRATCH/${SLURM_ARRAY_TASK_ID}.tmp $SCRATCH/${SLURM_ARRAY_TASK_ID}.tmp1 | sort | uniq | sed -e "s/\\_/\\:/g" > ${SCRATCH}/genebased/chr${SLURM_ARRAY_TASK_ID}_exclude_sets_${MAF}
done

