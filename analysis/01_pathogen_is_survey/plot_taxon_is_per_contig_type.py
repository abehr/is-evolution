import polars as pl
import plotly.express as px
from efm.config import data
from efm import dataset_definitions as defs, plot_utils as plt

samples = pl.scan_csv(data.ncbi.combined_final_samples)
isescan = pl.scan_csv(data.isescan.per_contig)
genomad = pl.scan_csv(data.genomad)

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

output_dir = data.output / '01_is_families_per_eskape'
output_dir.mkdir(exist_ok=True)
for tx in defs.eskape_order:
	# Filter to ESKAPE taxon
	pdf = df.filter(eskape=tx)

	if tx in taxa_with_untrustworthy_contig_classification:
		pdf = pdf.with_columns(contig_type=pl.lit('chromosome'))

	# Filter only to chromosome & plasmid contigs
	pdf = pdf.filter(pl.col('contig_type').is_in(['chromosome', 'plasmid']))

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
	contig_types = ['chromosome'] if tx in taxa_with_untrustworthy_contig_classification else ['chromosome', 'plasmid']
	contig_types = pl.DataFrame(dict(contig_type=contig_types))
	all_combos = all_samples.join(all_fams, how='cross').join(contig_types, how='cross')

	counts = (
		counts
		# Join with all possible combos & fill missing with 0
		.join(all_combos, on=col_group, how='right')
		.fill_null(0)
		.group_by('contig_type', 'display_fam')
		.agg(
			nIS_avg=pl.col('nIS').mean(),
			nIS_std=pl.col('nIS').std()
		)
		# Prevent lower error bar from being negative in the plot
		.with_columns(error_y_minus=pl
				.when(pl.col('nIS_std') > pl.col('nIS_avg')) # mean-std is negative
				.then('nIS_avg')
				.otherwise('nIS_std')
		)
	)

	# Sort fams by overall abundance, with 'other' always at the end
	fam_order = common_fams.to_list()
	if 'other' in counts['display_fam']: fam_order.append('other')

	# sort contigs & set their colors
	contig_order = ['chromosome', 'plasmid']
	contig_color = ['#DDA370', '#CA7328']

	fig = px.bar(
		counts, x='display_fam', y='nIS_avg', color='contig_type',
		error_y='nIS_std', error_y_minus='error_y_minus',
		category_orders=dict(display_fam=fam_order, contig_type=contig_order),
		color_discrete_sequence=contig_color,
		barmode='group',
		template='plotly_white'
	)

	plt.common_figure_style(fig, publish=True)
	fig.update_traces(error_y=dict(thickness=1))
	fig.update_yaxes(range=[0, 51], zerolinecolor='black', zerolinewidth=1)
	fig.update_layout(showlegend=False, yaxis_title=f'Mean copy per {tx} genome')
	fig.write_image(output_dir / f'{tx}.svg')
