import polars as pl
import plotly.express as px
from efm.config import data
from efm import plot_utils as plt

# Use the final set of filtered, quality-controlled public genomes with ISEScan results
samples = pl.read_csv(data.ncbi.combined_final_samples)

# Count ISL3 elements per genome
isl3 = (
	pl
	.scan_csv(data.isescan.per_contig)
	.filter(family='ISL3')
	.group_by('sample')
	.agg(nISL3=pl.col('count').sum())
	.collect()
)

# Create the dataframe to plot (ISL3 count per species)
df = (
	samples
	.select('sample', 'genus', 'species', 'genome_length')
	.join(isl3, on='sample', how='left')
	# Filter to the gram-positive cocci
	.filter(pl.col('genus').is_in(['Enterococcus', 'Staphylococcus', 'Streptococcus']))
	# Make sure samples with no ISL3 are still counted
	.fill_null(0, matches_supertype=True)
	.with_columns(nISL3_Mb=pl.col('nISL3') / (pl.col('genome_length')/1e6))
	
	# Also add a species display column
	.with_columns(sp=pl.format('{}. {}', pl.col('genus').str.slice(0, 1), 'species'))
)

# PLOT

sp_order = (
	df
	.group_by('genus', 'sp')
	.agg(pl.len().alias('count'), pl.col('nISL3_Mb').median()) # use median (vs mean) for parity with how boxplot render
	.sort(by=['genus', 'nISL3_Mb'], descending=[False, True])
)

fig = px.box(
	df, x='sp', y='nISL3_Mb', color='genus',
	color_discrete_map=plt.taxon_colors,
	category_orders=dict(sp=sp_order['sp'].to_list()),
	template='plotly_white'
)
# fig.update_traces(width=0.8, jitter=0.3, marker=dict(size=3))
plt.common_box_style(fig)
fig.update_traces(width=0.8)


#, line=dict(color='black', width=1), fillcolor='lightblue')
fig.update_layout(showlegend=False)
axis = dict(
	ticks='outside',
	ticklen=6,
	tickwidth=1,
	linecolor='black',
	color='black',
	showgrid=False,
	showline=True,
	zeroline=False
)

fig.update_yaxes(
	title=dict(text='ISL3 elements per Mb', font={'family':'Helvetica', 'size':12}),
	tickfont=dict(family='Helvetica', size=10),
	range=[-3, 33],
	**axis)
fig.update_xaxes(title='', tickangle=-45, tickfont=dict(family='Helvetica', style='italic', size=10), **axis)
plt.common_figure_style(fig, publish=True)


# Add number of genomes for each species
for _, sp, count, _ in sp_order.iter_rows():
	fig.add_annotation(
		x=sp, y=0, xref='x', yref='paper', text=str(count), 
		showarrow=False, font=dict(size=7)
	)

fig.update_layout(width=1000, height=400) # ensure that the tick labels render how we want before printing to file
fig.write_image(data.output / '02_ISL3_cocci.svg')