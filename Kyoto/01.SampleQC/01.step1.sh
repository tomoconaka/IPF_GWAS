#!/bin/bash

DATADIR=/home/nagasaki/workspace/pipeline/dist/v3/
SCRATCH=/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/
cat $DATADIR/output/verifyBamID2/verifyBamID2.summary.tsv | awk '$7 >= 0.03{print $1,$7}' | tail -n+2 > $SCRATCH/fail.verifyBamID2.sample
awk '$2 == "PAIR" && ($8 < 0.97 || $24 > 0.08) {print $1, $8, $24}' $DATADIR/gather.alignment_summary_metrics.tsv > $SCRATCH/fail.mappingrate0.97.chimera0.08.sample
cat $DATADIR/chr.idxstats.count.txt | tail -n+2 | awk '$2/$4 > 0.035 && $3/$4 >0.005' > $SCRATCH/fail.XXY.sample
cut -f1 /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/data/AGP3000sample.txt | tail -n+2 | sort | uniq > $SCRATCH/tmp1
awk '{print $1}' $SCRATCH/fail.verifyBamID2.sample $SCRATCH/fail.mappingrate0.97.chimera0.08.sample $SCRATCH/fail.XXY.sample | sort | uniq > $SCRATCH/tmp2
comm -23 $SCRATCH/tmp1 $SCRATCH/tmp2 > $SCRATCH/step1.tokeep.AGP3000.sample

cut -f1 /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/data/IPsample20220601.list | tail -n+2 | sort | uniq > $SCRATCH/tmp3
comm -23 $SCRATCH/tmp3 $SCRATCH/tmp2 > $SCRATCH/step1.tokeep.IP.sample

cat $SCRATCH/step1.tokeep.IP.sample | sed -e "s/PFKT2312/PFKT2312_wd/g" | sed -e "s/PFKT2334/PFKT2334_wd/g" | sed -e "s/PFKT2377/PFKT2377_wd/g"  > $SCRATCH/step1.tokeep.IP.rev.sample
