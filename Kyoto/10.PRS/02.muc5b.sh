SCRATCH_DIR=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/


plink2 \
  --pfile ${SCRATCH_DIR}/PRS/IPFPRS \
  --snp chr11:1219991_G/T \
  --export A \
  --out ${SCRATCH_DIR}/PRS/IPFMUC5B
