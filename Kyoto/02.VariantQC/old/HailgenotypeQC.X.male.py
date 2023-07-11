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

tmp = hl.import_vcf(str(args.arg1)+'.3.male.vcf.gz', array_elements_required=False, min_partitions=4, reference_genome='GRCh38', force_bgz=True)

mt = hl.variant_qc(tmp)

mt = mt.annotate_entries(AB = (mt.AD[1] / hl.sum(mt.AD) ))

mt = mt.filter_entries( (mt.GQ>=20) &
                 (mt.DP >= 5) )
dataset1 = mt.transmute_entries(GT = hl.if_else(mt.GT == hl.call(1), hl.call(1,1), hl.call(0,0)))

hl.export_vcf(dataset1, str(args.arg1)+'.4.male.vcf.bgz')
