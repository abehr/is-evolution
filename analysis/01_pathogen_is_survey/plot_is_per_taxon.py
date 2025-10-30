import polars as pl
import plotly.express as px
from efm.config import data
from efm import plot_utils as plt

samples = pl.scan_csv(data.ncbi.combined_final_samples)
isescan = pl.scan_csv(data.isescan.summary)

df = (
	samples
	.join(isescan, on='sample', how='left')
	.rename(dict(count='nIS'))
	.with_columns(nIS=pl.col('nIS').fill_null(0))
	.collect()
)

# Calculate IS per genome-Mb
df = df.with_columns(nIS_Mb=pl.col('nIS') / (pl.col('genome_length')/1e6))


# n = df['genus'].n_unique()
output_dir = data.output / '01_is_total_per_taxon'
output_dir.mkdir(exist_ok=True)
for genus in df['genus'].unique():
	pdf = df.filter(genus=genus)
	sp_order = pdf['species'].unique().sort().to_list()
	fig = px.box(
		pdf, x='species', y='nIS_Mb',
		category_orders=dict(species=sp_order),
		template='plotly_white'
	)
	fig.update_traces(
		fillcolor='#D5E9B8',
		line=dict(color='black', width=1),
		marker=dict(color='black', size=2)
	)
	fig.update_xaxes(tickangle=-45)
	fig.update_yaxes(range=[-8, 120], zeroline=False)# zerolinecolor='black', zerolinewidth=1)
	fig.update_layout(height=400)
	plt.common_figure_style(fig, publish=True)

	# We also want to include the number of genomes for each of these taxa:
	for sp, count in pdf['species'].value_counts().iter_rows():
		fig.add_annotation(
			x=sp, y=0, xref='x', yref='paper', text=str(count), 
			showarrow=False, font=dict(size=7)
		)

	fig.write_image(output_dir / f'{genus}.svg')