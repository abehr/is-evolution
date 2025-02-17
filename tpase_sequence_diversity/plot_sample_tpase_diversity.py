import polars as pl
import plotly.express as px
import sys

# Import plotting utils to add custom color scheme
sys.path.append('../..')
sys.path.append('../../eskape_tpase_clusters')
import display_clusters as plt
from util.efmpl import extract_name_from_product_col

def make_bar_plot(pdf, x_value, fam_order, text_col=None):
	fig = px.bar(pdf, x=x_value, y='fam', orientation='h', 
		category_orders=dict(fam=fam_order), 
		hover_data='product', 
		# text='pdt',
		template='plotly_white',
		color_discrete_sequence=['#779ECB'])
	
	
	if text_col:
		fig.update_traces(text=pdf[text_col])
		fig.update_traces(textfont=dict(color='white', size=14, family="Avenir"))
		fig.update_layout(uniformtext_minsize=14, uniformtext_mode='hide')

	plt.common_figure_style(fig, publish=True)
	fig.update_traces(marker_line_color="black", marker_line_width=0.1)
	return fig


# Inspect what the most common IS elements are, by nt or aa seqid
def avg_tpase_counts(tpase_df, col, avg_genome_length, fam_subset=None):
	avg = tpase_df.group_by(col).agg(
		fam=pl.col('fam').mode().first(),
		product=pl.col('product').mode().first(),
		count=pl.len(),
		avg_count=pl.len()/tpase_df.select('sample').n_unique(),
		avg_count_per_mb = (pl.len() / df.select('sample').n_unique()) / (avg_genome_length / 1e6)
	)
	if fam_subset:
		avg = avg.filter(pl.col('fam').is_in(fam_subset))

	# Note that the df must be sorted by count descending, in order for the bars to also be sorted by size.
	return avg.sort(by='avg_count', descending=True)


# ===== E. faecium exploration =====

df = pl.read_csv('data/efm/tpase_per_sample.csv', null_values='')
data = pl.read_csv('data/efm/samples.csv')
avg_genome_length = data.select('genome_size').to_series().mean()

# Mean tpase count / genome -> 113
df.group_by('sample').len().mean()

# Mean tpase count / genome per IS family -> tracks very well with what we expect from ISEScan actually.
record_counts_per_fam = df.group_by(['sample','fam']).len().group_by('fam').mean().sort(by='len', descending=True)
fam_order_1 = record_counts_per_fam.select('fam').to_series().to_list()
fam_order_2 = ['ISL3', 'IS256', 'IS30', 'IS6', 'IS3', 'IS110'] # order by count of most abundant tpase approximately

avg = avg_tpase_counts(df, 'nt_seqid', avg_genome_length)

avg = avg.with_columns(pdt=extract_name_from_product_col())

make_bar_plot(avg, 'avg_count', fam_order_2, 'pdt') # or avg count per mb