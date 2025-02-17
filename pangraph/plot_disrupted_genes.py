import polars as pl
import plotly.express as px
from ..figures import plot_utils as plt

other_fams = ['IS30','IS256','IS110','IS3','IS6']
info = f'ISL3 vs other commonn IS families ({', '.join(other_fams)})'

# ISL3: 11/100 intragenic
# Other: 20/51 intragenic
genes_disrupted = pl.DataFrame([
	dict(fam='ISL3', count=11, type='intra'),
	dict(fam='ISL3', count=89, type='inter'),
	dict(fam='Other', count=20, type='intra'),
	dict(fam='Other', count=31, type='inter')
])

fig = px.bar(
	genes_disrupted, x='count', y='fam', color='type',
	category_orders=dict(fam=['ISL3', 'Other'], type=['intra', 'inter']),
	color_discrete_sequence=['#E04B4B', '#6094C3'],
	template='plotly_white',
	title='Intragenic insertions (genes disrupted)'
	)

plt.common_figure_style(fig, publish=True)
fig.update_layout(title_pad_l=50, title_subtitle_text=info, showlegend=False)
fig.show()


# ================================================================================

ise = pl.concat([
	pl.read_ndjson('test_data/shc_st117/ISL3_intergenic_distance.jsonl').with_columns(fam=pl.lit('ISL3')),
	pl.read_ndjson('test_data/shc_st117/ISE_intergenic_distance.jsonl').with_columns(fam=pl.lit('Other'))
])

ise = ise.select('junction', 'fam', 'ise_dist', 'd_strand_parity')



fig = px.strip(
	ise, x='ise_dist', y='fam', color='d_strand_parity',
	category_orders=dict(fam=['ISL3', 'Other'], d_strand_parity=[True, False]),
	color_discrete_sequence=['#63BC6A', '#E04B4B'],
	template='plotly_white',
	title='Distance to nearest gene from IS insertion (from tpase orf 3\' end)'
	)

info = 'Distance to & relative orientation of nearest orf CDS for intergenic insertions. Same strand (green) or opposite strand (red).'
plt.common_figure_style(fig, publish=True)
fig.update_traces(jitter=0.8)
fig.update_layout(title_pad_l=50, title_subtitle_text=info, showlegend=False)