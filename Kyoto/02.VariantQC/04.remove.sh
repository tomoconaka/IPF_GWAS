SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch

rm ${SCRATCH}/${1}.vcf.gz
rm ${SCRATCH}/${1}.1.vcf.gz
rm ${SCRATCH}/${1}.{3..7}*.vcf*
rm -r ${SCRATCH}/${1}.final
