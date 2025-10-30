import polars as pl
import plotly.express as px

from efm.config import data
from efm import plot_utils as plt

samples = (
	pl
	.read_csv(data.shc.lr_meta)
	.filter(pl.col('genome').is_not_null()) # has Efm chromosome
	.select('patient', 'sample', 'barcode')
)

###############################
# Process "self" Sniffles
# (within-sample heterogeneity)
# and make a simple plot of this (Fig. 5B)
###############################

snf = pl.read_csv(data.lr_meta_sv.within_sample, separator='\t')

snf['svtype'].value_counts(sort=True)

pdf = (
	snf
	.group_by('sample').len()
	.join(samples.select('sample'), how='right', on='sample')
	.fill_null(0)
)
sample_order = pdf['sample'].sort().to_list()
fig = px.bar(
	pdf, x='sample', y='len',
	category_orders=dict(sample=sample_order),
	template='plotly_white'
)
plt.common_figure_style(fig, publish=True)
fig.write_image(data.output / '05_sv_within_sample.svg')

