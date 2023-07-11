#!/bin/bash
cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/Telseq

echo "ID telomere" > Telseq.summary
for i in $(ls *.result | sed -e "s/.result//g")
do
length=$(cat ${i}.result | sed -n 4p | cut -f7)
echo "${i} ${length}" >> Telseq.summary
done
bgzip Telseq.summary
