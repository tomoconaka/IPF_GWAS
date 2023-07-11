#!/bin/bash

GWASDIR=../../../../data/GWAS/
METADIR=../../../../Meta_scratch/
UKBDIR=../../../../UKB_scratch/

##Allen
zcat ${GWASDIR}/Allen_et_al/meta_gwas_5way_summary_stats.txt.gz | tail -n+2 | awk '{OFS="\t" ; print "chr"$1,int($2),int($2+1),$3}' > ${METADIR}/Allen_IPF.b37.bed
liftOver ${METADIR}/Allen_IPF.b37.bed ~/gpfs1/data/LiftOverchains/hg19ToHg38.over.chain ${METADIR}/Allen_IPF.b38.bed ${GWASDIR}/unMapped

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 N" > ${METADIR}/Allen_5way.txt
awk 'FNR==NR { m[$4]=$1; n[$4] = $2; next } \
($3 in m && $6 >= 0.001 && $6 <= 0.999) { OFS=" "; print m[$3]":"n[$3],$4,$5,$10,$11,$12,$6, $8}' ${METADIR}/Allen_IPF.b38.bed \
<(zcat ${GWASDIR}/Allen_et_al/meta_gwas_5way_summary_stats.txt.gz) |\
sed -e "s/chr//g" >> ${METADIR}/Allen_5way.txt
bgzip ${METADIR}/Allen_5way.txt

##FinnGen
echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 N" > ${METADIR}/FinnGen_IPF.txt
zcat ${GWASDIR}/finngen_R9_IPF.gz | tail -n+2 | awk -F "\t" '($11 >= 0.001 && $11 <= 0.999){OFS=" "; print $1":"$2,$3,$4,$9,$10,$7,$11,2018+373064}' >> ${METADIR}/FinnGen_IPF.txt
bgzip ${METADIR}/FinnGen_IPF.txt

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 N" > ${METADIR}/FinnGen_ILD.txt
zcat ${GWASDIR}/finngen_R7_ILD.gz | tail -n+2 | awk -F "\t" '($11 >= 0.001 && $11 <= 0.999){OFS=" "; print $1":"$2,$3,$4,$9,$10,$7,$11,3091+306063}' >> ${METADIR}/FinnGen_ILD.txt
bgzip ${METADIR}/FinnGen_ILD.txt

##UKB
zcat ${UKBDIR}/EUR_IPF.SAIGEformat.gz | tail -n+2 | awk '{OFS="\t"; print "chr"$1,$2,$2+1,$3}' | sed -e "s/chr23/chrX/g" | sed -e "s/chr24/chrX/g" > ${UKBDIR}/UKB_IPF.b37.bed
liftOver ${UKBDIR}/UKB_IPF.b37.bed ~/gpfs1/data/LiftOverchains/hg19ToHg38.over.chain ${UKBDIR}/UKB_IPF.b38.bed ${UKBDIR}/unMapped

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${METADIR}/UKB_IPF.txt
awk 'FNR==NR { m[$4]=$1; n[$4] = $2; next } \
($3 in m && $7 >= 0.3 && $6 >= 0.001 && $6 <= 0.999) { OFS=" "; print m[$3]":"n[$3],$4,$5,$10,$11,$14,$6,$7,$8}' ${UKBDIR}/UKB_IPF.b38.bed <(zcat ${UKBDIR}/EUR_IPF.SAIGEformat.gz) |\
sed -e "s/chrX/chr23/g" | sed -e "s/chr//g" >> ${METADIR}/UKB_IPF.txt
bgzip  ${METADIR}/UKB_IPF.txt

zcat ${UKBDIR}/EUR_ILD.SAIGEformat.gz | tail -n+2 | awk '{OFS="\t"; print "chr"$1,$2,$2+1,$3}' | sed -e "s/chr23/chrX/g" | sed -e "s/chr24/chrX/g" > ${UKBDIR}/UKB_ILD.b37.bed
liftOver ${UKBDIR}/UKB_ILD.b37.bed ~/gpfs1/data/LiftOverchains/hg19ToHg38.over.chain ${UKBDIR}/UKB_ILD.b38.bed ${UKBDIR}/unMapped

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${METADIR}/UKB_ILD.txt
awk 'FNR==NR { m[$4]=$1; n[$4] = $2; next } \
($3 in m && $7 >= 0.3 && $6 >= 0.001 && $6 <= 0.999) { OFS=" "; print m[$3]":"n[$3],$4,$5,$10,$11,$14,$6,$7,$8}' ${UKBDIR}/UKB_ILD.b38.bed <(zcat ${UKBDIR}/EUR_IPF.SAIGEformat.gz) |\
sed -e "s/chrX/chr23/g" | sed -e "s/chr//g" >> ${METADIR}/UKB_ILD.txt
bgzip  ${METADIR}/UKB_ILD.txt

#GBMI
echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${METADIR}/GBMI_AMR.txt
zcat ${GWASDIR}/IPF_Bothsex_amr_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$3,$4,$7,$8,$9,$6,1,$12+$13}' >> ${METADIR}/GBMI_AMR.txt
bgzip ${METADIR}/GBMI_AMR.txt

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${METADIR}/GBMI_AFR.txt
zcat ${GWASDIR}/IPF_Bothsex_afr_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$3,$4,$7,$8,$9,$6,1,$12+$13}' >> ${METADIR}/GBMI_AFR.txt
bgzip ${METADIR}/GBMI_AFR.txt


KyotoDIR=../../../../scratch/GWAS/
echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${KyotoDIR}/Kyoto_ILD.txt
cat ${KyotoDIR}/ILD.SAIGEformat.txt | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$4,$5,$9,$10,$15,$6,1,$7}' >> ${KyotoDIR}/Kyoto_ILD.txt

echo -e "SNP Allele1 Allele2 BETA SE p.value AF_Allele2 imputationInfo N" > ${KyotoDIR}/Kyoto_IPF.txt
cat ${KyotoDIR}/IPF.SAIGEformat.txt | tail -n+2 | awk '($6 >= 0.001 && $6 <= 0.999){print $1":"$2,$4,$5,$9,$10,$15,$6,1,$7}' >> ${KyotoDIR}/Kyoto_IPF.txt
