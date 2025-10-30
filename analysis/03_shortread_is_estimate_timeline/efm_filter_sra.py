import sys
import polars as pl
import plotly.express as px

# Custon filter accessions as specified in the Methods.

processed_entrez_results = sys.argv[1]
output_filtered_accessions_list = sys.argv[2]

missing_vals = ['missing', 'not applicable', 'not collected', 'not available' 'na', 'n/a', 'unknown', 'restricted access']

df = pl.read_csv(processed_entrez_results)

'''
# Possibly interesting: view filter/group by sample geolocation (though this isn't important for the current analysis)
df = (
	df
	.with_columns(country=pl.col('geo_loc_name').str.split(':').list.first())
	.with_columns(
		country=pl.when(
			pl.col('country').is_null() | 
			pl.col('country').str.to_lowercase().is_in(missing_vals))
		.then(pl.lit('MISSING'))
		.otherwise('country')
	)
)
px.histogram(df, x='year_disp', color='cov_disp', color_discrete_map=dict(low='red', OK='green', high='yellow'))
'''

# Also possibly interesting:
'''
Assigning host/env type
both can be null
when host is null, isolation_source often has info.

There are actually very very few samples where both host and isolation_source are both missing.
'''
missing = lambda col: (pl.col(col).is_null() | pl.col(col).str.to_lowercase().is_in(missing_vals))


# Instead:
# First, we'll def filter by paired & >20 est cov
include = (
	df
	.filter(library_layout='PAIRED')
	.filter(pl.col('est_cov') >= 20)
)

# Next, we'll filter to 500 results per year. 
def subsample_500(ydf):
	if len(ydf) <= 500: return ydf
	return ydf.sample(500, shuffle=True, seed=25)

include = pl.concat([subsample_500(x) for x in include.partition_by('year_disp')])


px.histogram(include, x='year_disp')


with open(output_filtered_accessions_list, 'w') as f:
	f.write('\n'.join(include['SRR']))