#!/bin/bash


cd /mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/genebased
~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNPID
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS ../../UKB_scratch/genebased/IPF.CADD20.0.001additive.SAIGEformat.txt.gz
PROCESS ../../Kyoto_scratch/genebased/IPF.CADD20.0.001additive.SAIGEformat.txt.gz

ANALYZE RANDOM
###
bgzip -c METAANALYSIS1.TBL > IPF.CADD20.0.001additive.b38.txt.gz


~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNPID
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS ../../UKB_scratch/genebased/IPF.CADD20.0.001dominant.SAIGEformat.txt.gz
PROCESS ../../Kyoto_scratch/genebased/IPF.CADD20.0.001dominant.SAIGEformat.txt.gz

ANALYZE RANDOM
###
bgzip -c METAANALYSIS1.TBL > IPF.CADD20.0.001dominant.b38.txt.gz

~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNPID
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS ../../UKB_scratch/genebased/IPF.CADD20.0.01recessive.SAIGEformat.txt.gz
PROCESS ../../Kyoto_scratch/genebased/IPF.CADD20.0.01recessive.SAIGEformat.txt.gz

ANALYZE RANDOM
###
bgzip -c METAANALYSIS1.TBL > IPF.CADD20.0.01recessive.b38.txt.gz

~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNPID
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS ../../UKB_scratch/genebased/IPF.LoF.0.01recessive.SAIGEformat.txt.gz
PROCESS ../../Kyoto_scratch/genebased/IPF.LoF.0.01recessive.SAIGEformat.txt.gz

ANALYZE RANDOM
###
bgzip -c METAANALYSIS1.TBL > IPF.LoF.0.01recessive.b38.txt.gz


~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNPID
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS ../../UKB_scratch/genebased/IPF.LoF.0.001dominant.SAIGEformat.txt.gz
PROCESS ../../Kyoto_scratch/genebased/IPF.LoF.0.001dominant.SAIGEformat.txt.gz

ANALYZE RANDOM
###
bgzip -c METAANALYSIS1.TBL > IPF.LoF.0.001dominant.b38.txt.gz

~/gpfs1/src/random-metal/bin/metal
###
SCHEME STDERR
AVERAGEFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER   SNPID
WEIGHT   N
ALLELE   Allele2 Allele1
EFFECT   BETA
FREQ     AF_Allele2
STDERR   SE
PVAL     p.value

PROCESS ../../UKB_scratch/genebased/IPF.LoF.0.001additive.SAIGEformat.txt.gz
PROCESS ../../Kyoto_scratch/genebased/IPF.LoF.0.001additive.SAIGEformat.txt.gz

ANALYZE RANDOM
###
bgzip -c METAANALYSIS1.TBL > IPF.LoF.0.001additive.b38.txt.gz


