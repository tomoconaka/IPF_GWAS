#!/bin/bash

GWASDIR=../../../../data/GWAS/
METADIR=../../../../Meta_scratch/
echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${METADIR}/GBMI_EAS.txt
zcat ${GWASDIR}/IPF_Bothsex_eas_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$3,$4,$7,$8,$9,$6,1,$12+$13}' >> ${METADIR}/GBMI_EAS.txt
bgzip ${METADIR}/GBMI_EAS.txt

zcat ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/GWASsummary_ILD_Japanese_SakaueKanai2020.auto.txt.gz | tail -n+2 | awk '{OFS="\t" ; print "chr"$2,int($3),int($3+1),$1}' > ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/BBJ_ILD.b37.bed
liftOver ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/BBJ_ILD.b37.bed ~/gpfs1/data/LiftOverchains/hg19ToHg38.over.chain ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/BBJ_ILD.b38.bed ${GWASDIR}/unMapped

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 N" > ${GWASDIR}/BBJ_ILD.txt
awk -F "\t" 'FNR==NR { m[$4]=$1; n[$4] = $2; next } \
($1 in m && $8 >= 0.001 && $8 <= 0.999) { OFS=" "; print m[$1]":"n[$1],$5,$6,$11,$12,$14,$8, $10}' ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/BBJ_ILD.b38.bed <(zcat ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/GWASsummary_ILD_Japanese_SakaueKanai2020.auto.txt.gz) |\
sed -e "s/chr//g" | sed -e "s/X:/23:/g" >> ${GWASDIR}/hum0197.v3.BBJ.ILD.v1/BBJ_ILD.txt

KyotoDIR=../../../../Kyoto_scratch/GWAS/
#echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${KyotoDIR}/Kyoto_ILD.txt
#cat ${KyotoDIR}/ILD.SAIGEformat.txt | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$4,$5,$9,$10,$15,$6,1,$7}' >> ${KyotoDIR}/Kyoto_ILD.txt

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${KyotoDIR}/Kyoto_IPF.txt
zcat ${KyotoDIR}/IPF.SAIGEformat.txt.gz | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$4,$5,$9,$10,$15,$6,1,$7}' >> ${KyotoDIR}/Kyoto_IPF.txt
bgzip ${KyotoDIR}/Kyoto_IPF.txt
