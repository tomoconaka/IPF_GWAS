SCRATCH_DIR=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/

for chr in {1..22}
do
tabix -p vcf ${SCRATCH_DIR}/${chr}.final.vcf.bgz &
done

for chr in {1..22}
do
bcftools view -Oz \
  -R ${SCRATCH_DIR}/IPF.extract \
  ${SCRATCH_DIR}/${chr}.final.vcf.bgz \
  -o ${SCRATCH_DIR}/PRS/IPFPRS.${chr}.vcf.gz
tabix -p vcf ${SCRATCH_DIR}/PRS/IPFPRS.${chr}.vcf.gz
done

ls ${SCRATCH_DIR}/PRS/IPFPRS.*.vcf.gz > ${SCRATCH_DIR}/PRS/IPFPRS.vcf.files

bcftools concat -Oz -f ${SCRATCH_DIR}/PRS/IPFPRS.vcf.files -a -o ${SCRATCH_DIR}/PRS/IPFPRS.vcf.gz

plink2 \
--vcf ${SCRATCH_DIR}/PRS/IPFPRS.vcf.gz \
--make-pfile \
--out ${SCRATCH_DIR}/PRS/IPFPRS

plink2 \
  --pfile ${SCRATCH_DIR}/PRS/IPFPRS \
  --score <(zcat ${SCRATCH_DIR}/../Meta_scratch/LeaveKyotoOut_meta_IPF.b38.txt.gz | awk 'NR != 1 {print "chr"$3"_"toupper($4)"/"toupper($5), toupper($4), $8}{print "chr"$3"_"toupper($5)"/"toupper($4), toupper($4), $8}') \
  --out ${SCRATCH_DIR}/PRS/IPFPRS
