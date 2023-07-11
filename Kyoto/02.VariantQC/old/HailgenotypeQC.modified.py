import hail as hl
from hail.plot import output_notebook, show
import bokeh
from pprint import pprint
import argparse
import os

os.chdir('/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/scratch/')

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

mt = mt.filter_rows((mt.variant_qc.gq_stats.mean > 10) & (mt.variant_qc.dp_stats.mean > 5) & (mt.variant_qc.call_rate > 0.8)  & (mtAll.variant_qc.p_value_hwe>5e-8))

hl.export_vcf(mt, str(args.arg1)+'.4.vcf.bgz')



