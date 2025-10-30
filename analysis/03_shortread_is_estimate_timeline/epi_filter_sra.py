import sys
import polars as pl
import plotly.express as px

# Custom filter accessions as specified in the Methods.

processed_entrez_results = sys.argv[1]
output_filtered_accessions_list = sys.argv[2]

df = pl.read_csv(processed_entrez_results)

# Note: could consider "nostril" or "nasal swab", but it's not crucial
skin = pl.col('envs').str.contains('skin') 

# Careful searching for words like "ear" which can be part of other words
skin2 = (
	pl.col('envs').is_in(['hand', 'foot', 'ear']) |
	pl.col('envs').str.contains_any(['nostril', 'antecubital', 'scalp', '_hand', '_foot', '_ear'])
)

hosp = pl.col('envs').str.contains_any(['blood', 'bacteremia', 'hospital'])

df = (
	df
	# null -> empty str, otherwise when you combine columns below with pl.format,
	# the entire thing becomes null if any one is null.
	.fill_null('')
	# concatenate all the relevant columns to determine source
	.with_columns(envs=pl.format(
		'{} {} {} {}', 
		'env_broad_scale', 'env_local_scale', 'env_medium', 'isolation_source'
	))
	# Strip whitespace & convert to lowercase
	.with_columns(envs=pl.col('envs').str.strip_chars().str.to_lowercase())
	.with_columns(
		# keep in mind that order matters here
		source=pl.when(skin).then(pl.lit('skin'))
				 .when(hosp).then(pl.lit('hosp'))
				 .when(skin2).then(pl.lit('skin2'))
				 .otherwise(pl.lit('other'))
	)
)

'''
Here, you can see that skin2 handles a lot of results from the same study in 2014.
After filtering skin & skin2, we're at a reasonable sample size. 
'''
px.histogram(df, x='year_disp', color='source')



# Filtering step:
include = (
	df
	.filter(pl.col('est_cov') >= 20)
	.filter(~pl.col('source').str.contains('skin'))
)

with open(output_filtered_accessions_list, 'w') as f:
	f.write('\n'.join(include['SRR']))