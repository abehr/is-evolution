import polars as pl
import plotly.express as px
from efm.config import data
from efm import dataset_definitions as defs, plot_utils as plt


samples = pl.read_csv(data.ncbi.combined_final_samples)
isescan = pl.read_csv(data.isescan.summary)

df = (
	samples
	.filter(pl.col('eskape').is_not_null())
	.join(isescan, on='sample', how='left')
	.rename(dict(count='nIS'))
	.with_columns(nIS=pl.col('nIS').fill_null(0)) # performative, basically
)

print(
	f'Public genomes from ESKAPEE genera with geNomad and ISEScan data:'
	f'\n{len(samples.filter(pl.col('collection').str.starts_with(defs.NCBI_NANOPORE)))} NCBI-longread'
	f'\n{len(samples.filter(pl.col('collection').str.starts_with(defs.NCBI_COMPLETE)))} NCBI-complete'
	f'\nOf these, {len(df)} are ESKAPEE organisms (shown in the plot).'
)

bar = (
	df # Calculate mean IS per genome-Mb
	.with_columns(nIS_Mb=pl.col('nIS') / (pl.col('genome_length')/1e6))
	.group_by('eskape') # for each ESKAPE taxon group
	.agg(
		mean=pl.mean('nIS_Mb'),
		std=pl.std('nIS_Mb'),
		nsamp=pl.len()
	)
	# You could use standard error as error bars (rather than stdev)
	# Casting doesn't seem necessary but w/e to be safe
	.with_columns(sem=pl.col('std')/pl.col('nsamp').cast(pl.Float64).sqrt())
)

fig = px.bar(bar, 
	x='mean', y='eskape', error_x='std',
	orientation='h', #text='nsamp',
	category_orders=dict(eskape=defs.eskape_order),
	template='plotly_white'
)

plt.common_figure_style(fig, y_italic=True, publish=True, small=False)
fig.update_traces(marker_color='#abd370', marker_line_color="#3a4726", marker_line_width=2)
fig.write_image(data.output / '01_eskape_IS.svg')