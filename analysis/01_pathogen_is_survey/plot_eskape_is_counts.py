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
	f'\n{len(samples.filter(pl.col("collection").str.starts_with(defs.NCBI_NANOPORE)))} NCBI-longread'
	f'\n{len(samples.filter(pl.col("collection").str.starts_with(defs.NCBI_COMPLETE)))} NCBI-complete'
	f'\nOf these, {len(df)} are ESKAPEE organisms (shown in the plot).'
)

# Calculate IS per genome-Mb
box = df.with_columns(nIS_Mb=pl.col('nIS') / (pl.col('genome_length')/1e6))

fig = px.box(
	box, x='nIS_Mb', y='eskape', orientation='h', 
	category_orders=dict(eskape=defs.eskape_order), 
	template='plotly_white'
)
plt.common_figure_style(fig, y_italic=True, publish=True, small=True)
plt.common_box_style(fig, optimize_smaller=False)
fig.update_traces(fillcolor='#D5E9B8', marker_color='#D5E9B8')

fig.write_image(data.output / '01_eskape_IS.svg')
box['sample', 'nIS_Mb', 'eskape'].write_excel(data.output / '01_eskape_IS.xlsx')