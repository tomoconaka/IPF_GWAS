#!/bin/bash

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/GWAS

cd ${SCRATCH}

echo "CHR POS ID Allele1 Allele2 AF_Allele2 N TEST BETA SE CHISQ LOG10P EXTRA" > ${1}.regenie
for i in {1..22} X PAR
do
tail -n+2 chr${i}_${1}.regenie >> ${1}.regenie
done
cat ${1}.regenie| awk '
NR==1{
sub("CHROM", "CHR");
sub("GENPOS", "POS");
sub("ID", "SNPID");
sub("ALLELE0", "Allele1");
sub("ALLELE1", "Allele2");
sub("A1FREQ", "AF_Allele2");
sub("INFO", "imputationInfo");
sub("N", "N")
for (i=1;i<=NF;i++) h[$i]=i; print $0,"SNP","p.value"
}
NR>1 && $h["LOG10P"]!="NA" {
print $0,$1":"$2,10^-$h["LOG10P"]
}' | bgzip -c > ${1}.SAIGEformat.txt.gz
