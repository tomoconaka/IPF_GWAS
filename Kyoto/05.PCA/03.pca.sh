#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
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

/home/tnakanishi/src/king -b ${SCRATCH}/grm_final.bed --unrelated --degree 2 
mv king* ${SCRATCH}

plink --bfile ${SCRATCH}/grm_final --keep ${SCRATCH}/kingunrelated.txt --make-bed --out ${SCRATCH}/grm_unrelated
plink --bfile ${SCRATCH}/grm_final --keep ${SCRATCH}/kingunrelated_toberemoved.txt --make-bed --out ${SCRATCH}/grm_related

flashpca --bfile ${SCRATCH}/grm_unrelated --outload ${SCRATCH}/grm_unrelated-loadings.pop \
--outmeansd ${SCRATCH}/grm_unrelated-meansd.pop --outpc ${SCRATCH}/grm_unrelated.pc --suffix .pop --ndim 20
flashpca --bfile ${SCRATCH}/grm_related --project --inmeansd ${SCRATCH}/grm_unrelated-meansd.pop --outproj ${SCRATCH}/grm_related.projections --inload ${SCRATCH}/grm_unrelated-loadings.pop -v

mv *.pop ${SCRATCH}
