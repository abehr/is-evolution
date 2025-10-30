import numpy as np, polars as pl
from scipy.stats import t as student_t
import plotly.express as px
import plotly.graph_objects as go
# note that statsmodels pkg is necessary for lowess smoothing, even though it is not explicitly called

from efm import plot_utils as plt

color_map = plt.IS_fam_color_map.copy()


def timeline_boxplot_binned(isc, fams, y_max=100):
	# px.box(isc, x='year_bin', y='est_copy', color='fam4', color_discrete_map=color_map, template='plotly_white')
	
	ybins = isc['year_bin'].unique().sort()
	
	box = go.Figure()
	for fam in fams:
		pdf = isc.filter(family=fam)
		for ybin in ybins:
			counts = pdf.filter(year_bin=ybin)['est_copy'].to_list()
			box.add_trace(go.Box(
				x=[fam]*len(counts),
				y=counts,
				marker=dict(color=color_map[fam], size=2),
				line=dict(color='black'),
				fillcolor=color_map[fam],
				opacity=0.7,
				boxpoints='all',
				offsetgroup=str(ybin),
				legendgroup=str(ybin),
				pointpos=-2, # move the points pretty far to the left
				jitter=0.5 # increase jitter so you can see them better
			))
	box.update_layout(
		boxmode='group', template='plotly_white', showlegend=False,
		boxgap=0.2, # keep a small gap between the family groups
		boxgroupgap=0.5 # increase the gap between the boxes within a family group to make room for the points
		)
	
	plt.common_figure_style(box, publish=True, small=False)

	# Y axis information
	outliers = isc.filter(pl.col('family').is_in(fams), pl.col('est_copy') > y_max)
	ytitle = f'Estimated IS copy per genome ({len(outliers)} outlier values >{y_max} not shown)'
	box.update_yaxes(range=[-5, y_max], title=ytitle)

	# X axis information
	bin_counts = (
		isc['sample', 'year_bin'].unique()
		.group_by('year_bin').len()
		.sort('year_bin')
		.with_columns(disp=pl.format('{}-{} ({} samples)', pl.col('year_bin')-4, 'year_bin', 'len'))
	)
	xtitle = 'Sample year bins: ' + ', '.join(bin_counts['disp'])

	box.update_xaxes(title=xtitle)
	return box



def mean_t_SE_CI(vals, conf=.95):
	n = vals.size
	if n < 2: return np.nan, np.nan

	mean = np.mean(vals)
	sd = np.std(vals, ddof=1)
	se = sd / np.sqrt(n)
	t_crit = student_t.ppf(0.5+conf/2.0, n-1)
	half = se*t_crit
	return mean, half


# This can just take the "isc" plot. We'll just implement mean for now
def timeline_lowess_sliding_window(df, fams, window=2, lowess_smoothing=.25, fam_labels=False):
	aggs = []
	for year in range(df['year'].min()+window, df['year'].max()-window+1):
		per_fam = df.filter(pl.col('year').is_between(year-window, year+window)).group_by('family').agg('est_copy')
		for fam in fams:
			vals = per_fam.filter(family=fam)['est_copy'].item().to_numpy()
			mean, mse = mean_t_SE_CI(vals)
			aggs.append(dict(year=year, family=fam, mean=mean, mse=mse))
	
	aggs = pl.DataFrame(aggs)

	fig = px.scatter(
		aggs, x='year', y='mean', error_y='mse',
		color='family', color_discrete_map=color_map,
		trendline='lowess', trendline_options=dict(frac=lowess_smoothing), 
		template='plotly_white'
	)

	# Decrease the opacity of the dots & error bars.
	for trace in fig.data:
		faded_color = plt.hex2rgb(color_map[trace.name], a=0.5)
		if trace.mode == 'markers': # points & error bars are marker mode
			trace.update(
				marker=dict(size=5), 
				error_y=dict(color=faded_color, thickness=1, width=3)
			)
		elif trace.mode == 'lines': # LOWESS curves are line mode
			trace.update(line=dict(width=3))
			
			# Add the name of the IS family to the right of the end of the line
			if fam_labels:
				plt.add_family_name(fig, trace.name, color_map[trace.name], trace.y[-1])

	plt.common_figure_style(fig, publish=True, small=False)
	
	# Hide the legend and add some spacing, if we're plotting the fam name labels
	if fam_labels: fig.update_layout(showlegend=False, margin=dict(r=150))
	fig.update_xaxes(range=[aggs['year'].min()-1, aggs['year'].max()+1])
	fig.update_yaxes(range=[-2, int(aggs['mean'].max()+2)])

	return fig