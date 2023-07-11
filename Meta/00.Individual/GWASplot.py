from scipy.stats import uniform
from scipy.stats import randint
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

def manhattan(INPUT, OUTDIR, PHENO, PVAL):
	tmp = pd.read_table(INPUT, sep=" ")
	df = tmp.dropna(subset=['LOG10P'])
	df['ind'] = range(len(df))
	df_grouped = df.groupby(('CHR'))
	fig = plt.figure(figsize=(10,5))
	fig.suptitle(PHENO, size=18)
	ax = fig.add_subplot(111)
	colors = ['gray','black']
	x_labels = []
	x_labels_pos = []
	for num, (name, group) in enumerate(df_grouped):
		group.plot(kind='scatter', x='ind', y='LOG10P',color=colors[num % len(colors)], ax=ax)
		x_labels.append(name)
		x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels, rotation=45)
	ax.set_xlim([0, max(len(df), 10)])
	ax.set_ylim([0, max(max(df.LOG10P)+2, 9)])
	ax.set_xlabel('Chromosome')
	ax.axhline(-np.log10(PVAL), ls = "-.", color = "navy")
	fig.savefig(OUTDIR + PHENO + '.png')

def qqplot(INPUT, OUTDIR, PHENO):
	import hail as hl
	from hail.plot import output_notebook, show
	import bokeh
	from bokeh.io import export_png
	from pprint import pprint
	tmp = pd.read_table(INPUT, sep=" ")
	#df = tmp.dropna(subset=['LOG10P'])
	df = tmp.loc[:,['p.value']]
	a = hl.Table.from_pandas(df)
	p = hl.plot.qq(a['p.value'], title=PHENO)
	export_png(p, filename= OUTDIR + PHENO + '.qqplot.png')


def main(args):
	if args.manhattan:
		manhattan(args.input, args.outdir, args.pheno, args.pval)	

	if args.qqplot:
		qqplot(args.input, args.outdir, args.pheno)
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', type=str, required=True)
	parser.add_argument('--outdir', type=str, required=True)
	parser.add_argument('--pheno', type=str, required=True)
	parser.add_argument('--pval', type=float, required=False)
	parser.add_argument('--manhattan', action='store_true')
	parser.add_argument('--qqplot', action='store_true')
	args = parser.parse_args()
	main(args)
