#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=5
#SBATCH --job-name=vep
#SBATCH --mail-user=tnakanishi@genome.med.kyoto-u.ac.jp
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -p ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
##SBATCH --array=2

DATADIR=/home/nagasaki/workspace/pipeline/output/VQSR.v3
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

vep38="singularity run --bind /mnt /mnt/file1/pub/data/vep/ensembl-vep_release_104.3b.sif vep \
--offline --dir_cache /mnt/gpfs1/home/tkawa/kgwast/work/all/20210407/vep-docker --species homo_sapiens \
--assembly GRCh38 --force_overwrite --tab --compress_output bgzip --everything "

LOFTEE_PATH0=/mnt/file1/pub/data/vep/plugin/loftee
ANCESTOR_PATH=$LOFTEE_PATH0/GRCh38/human_ancestor.fa.gz
GERP_PATH=$LOFTEE_PATH0/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
CONSERVATION_PATH=$LOFTEE_PATH0/GRCh38/loftee.sql
REF_FASTA=/mnt/file1/pub/data/vep/homo_sapiens/104_GRCh38/hg38.fa.gz

## singularityの中
LOFTEE_PATH1=/opt/vep/Plugins/loftee

SLURM_ARRAY_TASK_ID=X
## 実行する
$vep38 -o $SCRATCH/LoF.chr${SLURM_ARRAY_TASK_ID}.txt.gz -i ${SCRATCH}/${SLURM_ARRAY_TASK_ID}.2.vcf.gz \
--fa $REF_FASTA \
--dir_plugins $LOFTEE_PATH1 \
--plugin LoF,loftee_path:$LOFTEE_PATH1,human_ancestor_fa:$ANCESTOR_PATH,conservation_file:$CONSERVATION_PATH,gerp_bigwig:$GERP_PATH

