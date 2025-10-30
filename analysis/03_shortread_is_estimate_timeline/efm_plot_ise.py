import polars as pl
from ise_plots import timeline_boxplot_binned, timeline_lowess_sliding_window
from efm.config import data

# Short-read IS count estimates produced by count_ise script
efm_is_est_counts_fp = data.timeline / 'efaecium/est_is_counts.csv'
if not efm_is_est_counts_fp.is_file():
	raise FileNotFoundError(
		'E. faecium estimated IS counts data does not exist '
		'-- run count_ise.py to generate it.'
	)

isc = (
	pl
	.read_csv(efm_is_est_counts_fp)
	.rename(dict(year_disp='year'))
	.filter(pl.col('year') != 2025) # ignore samples from 2025 (partial year)
	.with_columns(
		# Group together all pre-1995 samples bc there are so few of them. 
		year = pl.when(pl.col('year') < 1995).then(1994).otherwise('year'),
		# Add 5-year bins
		year_bin = pl.when(pl.col('year') <= 1995).then(1995).otherwise(((pl.col('year') - 1996) // 5) * 5 + 2000)
	)
)

efm_common_fams = ['ISL3', 'IS30', 'IS256', 'IS3', 'IS200/IS605', 'IS110']

### 5-year bins boxplot
box = timeline_boxplot_binned(isc, efm_common_fams)
box.write_image(data.output / '03_efm_timeline_boxplot.svg')

### Sliding window timeline with 95% CI
line = timeline_lowess_sliding_window(isc, efm_common_fams, window=2, lowess_smoothing=.25)
line.write_image(data.output / '03_efm_timeline_trendline.svg')
