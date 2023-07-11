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

tmp = hl.import_vcf(str(args.arg1)+'.3.vcf.gz', array_elements_required=False, min_partitions=4, reference_genome='GRCh38', force_bgz=True)

mt = hl.variant_qc(tmp)

mt = mt.annotate_entries(AB = (mt.AD[1] / hl.sum(mt.AD) ))

mt = mt.filter_entries( (mt.GQ>=20) &
                 (mt.DP >= 10) &
                 ((mt.GT.is_hom_ref() & (mt.AB <= 0.1)) |
                        (mt.GT.is_het() & (mt.AB >= 0.25) & (mt.AB <= 0.75)) |
                        (mt.GT.is_hom_var() & (mt.AB >= 0.9)))) 

mt = mt.filter_rows(mt.variant_qc.p_value_hwe>1e-6)
LCR_intervals = hl.import_locus_intervals('/home/tnakanishi/gpfs1/data/1000G/LCR/b38.bed', reference_genome='GRCh38')
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.filter_rows(mt.in_LCR, keep=False)

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



