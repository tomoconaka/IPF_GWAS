
#bash 06.convert_regenie2SAIGE.sh ILD CADD20.0.001 dominant

SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/genebased/

cd ${SCRATCH}

echo "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA" > ${1}.${2}${3}.regenie
for i in {1..22} X PAR
do
tail -n+2 chr${i}_${3}_${1}.regenie | grep ${2} >> ${1}.${2}${3}.regenie
done
cat ${1}.${2}${3}.regenie| awk '
NR==1{
sub("CHROM", "CHR");
sub("GENPOS", "POS");
sub("ID", "SNPID");
sub("ALLELE0", "Allele1");
sub("ALLELE1", "Allele2");
sub("A1FREQ", "AF_Allele2");
sub("N", "N")
for (i=1;i<=NF;i++) h[$i]=i; print $0,"SNP","p.value"
}
NR>1 && $h["LOG10P"]!="NA" {
print $0,$1":"$2,10^-$h["LOG10P"]
}' | bgzip -c > ${1}.${2}${3}.SAIGEformat.txt.gz
rm ${1}.${2}${3}.regenie
