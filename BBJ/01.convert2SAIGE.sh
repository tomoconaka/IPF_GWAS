
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/BBJ_scratch/

cd ${SCRATCH}

echo "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA" > ILD.regenie
for i in {1..22} X
do
for j in p q
do
if test -f "chr${i}${j}_ILD.regenie.gz"; then
    zcat chr${i}${j}_ILD.regenie.gz | tail -n+2 >> ILD.regenie
fi
done
done

cat ILD.regenie| awk '
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
}' | bgzip -c > ILD.SAIGEformat.txt.gz



