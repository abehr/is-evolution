import numpy as np
import polars as pl
import plotly.express as px
import xlsxwriter
from efm.config import data
from efm import dataset_definitions as defs, plot_utils as plt

samples = pl.scan_csv(data.ncbi.combined_final_samples)
isescan = pl.scan_csv(data.isescan.per_contig)
genomad = pl.scan_csv(data.genomad)

# Include E. faecalis in this analysis, even though it is no longer considered ESKAPE(E)
samples = samples.with_columns(
	eskape=pl
		.when(species='faecalis')
		.then(pl.lit('E. faecalis'))
		.otherwise('eskape')
)

df = (
	samples
	.filter(pl.col('eskape').is_not_null())
	.join(genomad, on='sample', how='left')
	.join(isescan, on=['sample', 'contig'], how='left')
	.fill_null(strategy='zero', matches_supertype=True)
	.rename(dict(count='nIS'))
	.collect()
)

# contig_types = pl.DataFrame(dict(contig_type=['chromosome', 'plasmid'])) # only consider these contig types
col_group = ['sample', 'display_fam', 'contig_type']

# Do not factor in contig type for those taxa where geNomad result is often unconfident at the contig level
taxa_with_untrustworthy_contig_classification = ['E. coli', 'Enterobacter spp.', 'K. pneumoniae']

source_data = []
output_dir = data.output / '01_is_families_per_eskape'
output_dir.mkdir(exist_ok=True)
for tx in defs.eskape_order + ['E. faecalis']:
	# Filter to ESKAPE taxon
	pdf = df.filter(eskape=tx)

	# sort contigs & set their colors
	contig_order = ['chromosome', 'plasmid']
	# contig_color = ['#DDA370', '#CA7328']
	contig_color = ['#A97353', '#6E1800']

	if tx in taxa_with_untrustworthy_contig_classification:
		pdf = pdf.with_columns(contig_type=pl.lit('all contigs'))
		contig_order = ['all contigs']
		contig_color = ['#F4723A']

	# Filter only to chromosome & plasmid contigs
	pdf = pdf.filter(pl.col('contig_type').is_in(['all contigs', 'chromosome', 'plasmid']))

	# Get the 10 most common IS families and set others to 'other'
	common_fams = pdf.group_by('family').agg(pl.sum('nIS')).sort(by='nIS', descending=True)
	common_fams = common_fams['family'][:10]
	pdf = pdf.with_columns(
		display_fam=pl
			.when(pl.col('family').is_in(common_fams.implode()))
			.then('family')
			.otherwise(pl.lit('other'))
	)

	# IS counts per sample (broken down by fam & contig type)
	counts = pdf.group_by(col_group).agg(pl.col('nIS').sum())

	# However, these counts omit zeros, which would artificially inflate the
	# averages. Therefore, join with all possible combos & fill missing with 0.
	all_samples = pdf.select('sample').unique()
	all_fams = pdf.select('display_fam').unique()
	contig_types = ['all contigs'] if tx in taxa_with_untrustworthy_contig_classification else ['chromosome', 'plasmid']
	contig_types = pl.DataFrame(dict(contig_type=contig_types))
	all_combos = all_samples.join(all_fams, how='cross').join(contig_types, how='cross')

	counts = (
		counts
		# Join with all possible combos & fill missing with 0
		.join(all_combos, on=col_group, how='right')
		.fill_null(0)
	)

	# Sort fams by overall abundance, with 'other' always at the end
	fam_order = common_fams.to_list()
	if 'other' in counts['display_fam']: fam_order.append('other')

	fig = px.box(
		counts, x='display_fam', y='nIS', color='contig_type', points='outliers',
		category_orders=dict(display_fam=fam_order, contig_type=contig_order),
		color_discrete_sequence=contig_color,
		template='plotly_white'
	)
	plt.common_figure_style(fig, publish=True)
	plt.common_box_style(fig)

	# Bring the groups a bit closer together to decrease horizontal whitespace
	fig.update_layout(boxgroupgap=0.1, boxgap=0.2)

	# A bit hacky, but basically we want to pretend to plot the outliers (on the boxplot) 
	# so that the axes are sized correctly and so that the ticks carry the intended meaning 
	# (1.5*IQR rather than min/max). But then we don't want to show them, because we'll plot 
	# all points on a separate trace (strip plot):
	fig.update_traces(marker=dict(color='white'))
	points = px.strip(counts, x='display_fam', y='nIS', color='contig_type', stripmode='group')
	rng = np.random.default_rng(0)
	for trace in points.data:
		y_force_points = trace.y + rng.uniform(-1e-3, 1e-3, size=len(trace.y))
		trace.update(jitter=0.6, opacity=0.5, marker=dict(size=2, color='#333333'), y=y_force_points)

	fig.add_traces(points.data)

	# For the updated boxplot version you actually don't want to have black & 1-width zeroline
	# because it is visually indistinguishable from the box whiskers.
	fig.update_layout(showlegend=False, yaxis_title=f'Mean copy per {tx} genome')
	fig.write_image(output_dir / f'{tx}.svg')
	source_data.append(counts)


# Write out plot source data
with xlsxwriter.Workbook(output_dir / 'source_data.xlsx') as xl:
	for tx,df in zip(defs.eskape_order + ['E. faecalis'], source_data):
		df.write_excel(xl, xl.add_worksheet(tx))