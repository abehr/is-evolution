import sys
import polars as pl
from pathlib import Path


library_layout_rank = dict(PAIRED=1, SINGLE=0)
instrument_rank = ['NovaSeq', 'HiSeq', 'HiScanSQ', 'NextSeq', 'MiSeq', 'MiniSeq'] # Genome Analyzer would be last
instrument_rank = {x:i for i,x in enumerate(instrument_rank)}

def main():
	if len(sys.argv) != 5:
		sys.exit(
			'Usage: python process_entrez_results.py <approx_genome_size_mb> <samples_ndjson> <experiments_ndjson> <output_csv>\n'
			'E.g.: python process_entrez_results.py 2.8 efm_samples.jsonl efm_experiments.jsonl processed_entrez_results.csv'
		)

	est_genome_size = float(sys.argv[1])*1E6
	samples = pl.read_ndjson(sys.argv[2]).unique() # file would have same nrows as e, but some are dupes.
	experiments = pl.read_ndjson(sys.argv[3])
	output_csv = Path(sys.argv[4])
	assert not output_csv.exists(), f'Output CSV {output_csv} already exists!'
	print(
		f'{len(samples)} samples ({sys.argv[2]})\n'
		f'{len(experiments)} experiments ({sys.argv[3]})\n'
		f'Est Genome Size: {est_genome_size/1E6} Mb'
	)

	df = merge(samples, experiments)
	df = with_valid_collection_year(df)
	df = with_est_coverage(df, est_genome_size)

	## We can probably ignore any non-paired-end data, because there's not enough to warrant it. Note:
	df.filter(pl.col('library_layout') == 'SINGLE', pl.col('est_cov') > 20) # empty for Staph epi
	df.filter(pl.col('library_layout') == 'SINGLE', pl.col('est_cov') > 20, pl.col('year_disp') < 2005) # only three are pre-2005 for Efm

	## We can also probably ignore any data that's spread across multiple SRR accessions (due to size).
	## Almost all of these are too-high coverage, which we could downsample, but they're also mostly from newer years,
	## where there is already sufficient sample size.
	df.filter(pl.col('SRR').str.contains(','))['year'].value_counts()

	
	# For both Efm & Staph epi, it's actually quite rare to find too-low coverage pre-2005, but is common to find too-high coverage.
	# Thus we will want to rarify down.
	plot_coverage(df)

	df.write_csv(output_csv)



'''
Note that SRR can have multiple values (comma-separated), because there can be 
multiple runs for a given experiment. This can be due to filesize. You can just combine them 
before assembly. I think you also can from SRA.
'''	

# Filter to only one experiment per sample. There are sometimes multiple (e.g. several runs
# were done on the same sample, sometimes with different sequencing tech)
def merge(samples, experiments):
	assert set(samples.columns).intersection(set(experiments.columns)) == {'SRS'}, 'Overlapping columns (cannot merge)'
	return (
		experiments
		.with_columns(pl.col('instrument_model').str.replace('Genome Analyzer', 'GenomeAnalyzer', literal=True))
		.with_columns(split=pl.col('instrument_model').str.strip_prefix("Illumina ").str.split(' '))
		.with_columns(
			model_name=pl.col('split').list.first(),
			model_detail=(
				pl.col('split').list.get(1, null_on_oob=True)
				.replace('X', '1000000') # "X" series instruments categorized highest
				.cast(int, strict=False)
			)
		)
		.drop('split')
		.with_columns(
			layout_rank=pl.col('library_layout').replace_strict(library_layout_rank, default=0),
			instrument_rank=pl.col('model_name').replace_strict(instrument_rank, default=99)
		)
		.sort(
			by=[
				"SRS", 
				"layout_rank",     # PAIRED > SINGLE
				"instrument_rank", # NovaSeq, HiSeq, ..., unspecified
				"model_detail",    # Highest > Lowest (X = arbitrarily high number)
				"bases"            # Highest > Lowest
			],
			descending=[False, True, False, True, True]
		)
		.unique(subset='SRS', maintain_order=True, keep='first')
		.drop('layout_rank', 'instrument_rank') # drop sorting columns
		.join(samples, on='SRS', how='left') # Optional to join...
	)

'''
Collection date: there are occasionally ones with date ranges, e.g. '2012-01-01/2015-06-30'.
These are actually very common -- out of ~9.5k Staph epi samples, 3.5k are a range.
Partially bc the majority of samples are from only a few studies. 
Let's do it like this: If the date range is ~3 years or less, take the middle, because the resulting
year won't be off by >~1. Otherwise ignore, because the date window is too wide to get an accurate year.
'''
def with_valid_collection_year(df):
	return (
		df
		# Get start/end date. These will be the same if the collection date is not a range.
		.with_columns(ysplit=pl.col('collection_date').fill_null('').str.split('/'))
		.with_columns(
			ystart=pl.col('ysplit').list.first().str.split('-').list.first().cast(int, strict=False),
			yend=pl.col('ysplit').list.last().str.split('-').list.first().cast(int, strict=False)
		)
		.with_columns(ydiff=pl.col('yend') - pl.col('ystart'))
		.filter(pl.col('ydiff').is_between(0, 3))
		.with_columns(year=pl.col('ystart') + (pl.col('ydiff')/2).cast(int))
		.with_columns(year_disp=pl.when(pl.col('year') < 1990).then(pl.lit(1990)).otherwise('year'))
		.drop('ysplit', 'ystart', 'yend', 'ydiff') # get rid of intermediate columns
	)

# Estimate coverage and assign est threshold to get a visual sense of sample size and what we might cut
def with_est_coverage(df, est_genome_size):
	return (
		df
		.with_columns(est_cov=pl.col('bases') / est_genome_size)
		.with_columns(
			cov_disp=pl.when(pl.col('est_cov') < 20)
						.then(pl.lit('low'))
						.when(pl.col('est_cov') > 200)
						.then(pl.lit('high'))
						.otherwise(pl.lit('OK'))
		)
	)

def plot_coverage(df):
	import plotly.express as px
	cov_plt = px.histogram(df, x='year_disp', color='cov_disp', color_discrete_map=dict(low='red', OK='green', high='yellow'))
	#cov_plt.write_image('cov.png')
	return cov_plt


if __name__ == '__main__':
	main()