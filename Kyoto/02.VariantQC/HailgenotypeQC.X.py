import hail as hl
from hail.plot import output_notebook, show
import bokeh
from pprint import pprint
import argparse
import os

os.chdir('/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Kyoto_scratch/')

hl.init(local='local[40]',min_block_size=128, tmp_dir='/home/tnakanishi/gpfs1/tmp/')

parser = argparse.ArgumentParser()
parser.add_argument('arg1')

args = parser.parse_args()

mt = hl.import_vcf(str(args.arg1)+'.3.vcf.gz', array_elements_required=False, min_partitions=4, reference_genome='GRCh38', force_bgz=True)

male = hl.import_table('all.sample.male', no_header=True)
female = hl.import_table('all.sample.female', no_header=True)
male1 = male.f0.collect()
males_to_keep = set(male1)
female1 = female.f0.collect()
females_to_keep = set(female1)

LCR_intervals = hl.import_locus_intervals('/home/tnakanishi/gpfs1/data/1000G/LCR/b38.bed', reference_genome='GRCh38')
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.filter_rows(mt.in_LCR, keep=False)

mt_male = mt.filter_cols(hl.literal(males_to_keep).contains(mt.s))
mt_female = mt.filter_cols(hl.literal(females_to_keep).contains(mt.s))

mt_male = hl.variant_qc(mt_male)
mt_female = hl.variant_qc(mt_female)

mt_male = mt_male.annotate_entries(AB = (mt_male.AD[1] / hl.sum(mt_male.AD) ))
mt_male = mt_male.filter_entries( (mt_male.GQ>=20) &
                 (mt_male.DP >= 5))

mt_male = mt_male.transmute_entries(GT = hl.if_else(mt_male.GT == hl.call(1), hl.call(1,1), mt_male.GT))
mt_male = mt_male.transmute_entries(GT = hl.if_else(mt_male.GT == hl.call(0), hl.call(0,0), mt_male.GT))

mt_female = mt_female.annotate_entries(AB = (mt_female.AD[1] / hl.sum(mt_female.AD) ))

mt_female = mt_female.filter_entries( (mt_female.GQ>=20) &
                 (mt_female.DP >= 10) &
                 ((mt_female.GT.is_hom_ref() & (mt_female.AB <= 0.1)) |
                        (mt_female.GT.is_het() & (mt_female.AB >= 0.25) & (mt_female.AB <= 0.75)) |
                        (mt_female.GT.is_hom_var() & (mt_female.AB >= 0.9)))) 

mt_female = mt_female.filter_rows(mt_female.variant_qc.p_value_hwe>1e-6)

mt = mt_male.union_cols(mt_female)
hl.export_vcf(mt, str(args.arg1)+'_malefemale.vcf.bgz')

mt = hl.import_vcf(str(args.arg1)+'_malefemale.vcf.bgz', array_elements_required=False, min_partitions=4, reference_genome='GRCh38', force_bgz=True)
mt = hl.variant_qc(mt)
cases = hl.import_table('SampleQC.tokeep.ILD.sample', no_header=True)
controls = hl.import_table('SampleQC.tokeep.AGP3000.sample', no_header=True)
cases1 = cases.f0.collect()
cases_to_keep = set(cases1)
controls1 = controls.f0.collect()
controls_to_keep = set(controls1)

cases_mt = mt.filter_cols(hl.literal(cases_to_keep).contains(mt.s))
controls_mt = mt.filter_cols(hl.literal(controls_to_keep).contains(mt.s))

cases_mt = hl.variant_qc(cases_mt)
controls_mt = hl.variant_qc(controls_mt)
cases_mt = cases_mt.filter_rows(cases_mt.variant_qc.call_rate > 0.8)
controls_mt = controls_mt.filter_rows(controls_mt.variant_qc.call_rate > 0.8)

final_mt = cases_mt.union_cols(controls_mt)

hl.export_vcf(final_mt, str(args.arg1)+'.final.vcf.bgz')



