import polars as pl
from ise_plots import timeline_lowess_sliding_window
from efm.config import data

# Short-read IS count estimates produced by count_ise script
epi_is_est_counts_fp = data.timeline / 'staph_epi/est_is_counts.csv'
if not epi_is_est_counts_fp.is_file():
	raise FileNotFoundError(
		'S. epidermidis estimated IS counts data does not exist '
		'-- run count_ise.py to generate it.'
	)

isc = (
	pl
	.read_csv(epi_is_est_counts_fp)
	.rename(dict(year_disp='year'))
	.filter(pl.col('year') != 2025) # ignore samples from 2025 (partial year)
	.with_columns(
		# Group together all pre-2000 samples bc there are so few of them. 
		year = pl.when(pl.col('year') < 2000).then(1999).otherwise('year'),
		# Add 5-year bins
		year_bin = pl.when(pl.col('year') <= 2000).then(2000).otherwise(((pl.col('year') - 2001) // 5) * 5 + 2005)
	)
)


fams = isc.group_by('family').agg(pl.sum('count')).sort('count', descending=True)['family'].to_list()
fams.remove('other')
fams = fams + ['other']

fig = timeline_lowess_sliding_window(isc, fams, window=2, lowess_smoothing=.5)
fig.update_yaxes(range=[0,13])
fig.write_image(data.output / '03_epi_timeline_trendline.svg')