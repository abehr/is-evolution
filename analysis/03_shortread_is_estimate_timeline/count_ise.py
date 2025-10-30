import sys
import polars as pl
from efm.config import data

'''
Use ISEScan and assembly coverage info to estimate IS counts per genome from short-read data. 

This script can be run after completing the pre.nf and post.nf short-read analysis pipelines, 
concatenating the results, and performing QC. This script will filter based on a qc_pass column.
See our Zenodo project for example file formats if you want to run this on your own data. 
'''

def main():
	if len(sys.argv) != 2:
		sys.exit('count_ise.py <source_files_dir_name>')

	taxon_name = sys.argv[1]
	
	wdir = data.timeline / taxon_name
	if not wdir.is_dir():
		raise FileNotFoundError(f'Input directory {wdir} does not exist')

	contig_stats = pl.read_csv(wdir / 'assembly_contigs_parsed.csv')
	isescan_sum = pl.read_csv(wdir / 'isescan_sum.csv')
	sample_qc = pl.read_csv(wdir / 'assembly_qc.csv')

	# Validate that the input files have the information we need
	assert contig_stats.columns == ['sample', 'contig', 'length', 'cov']
	assert all(x in isescan_sum.columns for x in ('sample', 'contig', 'family', 'count'))
	assert all(x in sample_qc.columns for x in ('sample', 'year_disp', 'qc_pass'))
	
	# Output file is named after the input dir name
	output_csv = wdir / f'est_is_counts.csv'


	# "Final" list of samples = intersection of samples with ISEScan results 
	# & samples that passed assembly/species QC
	isescan_samples = isescan_sum.filter(family='total')['sample'].unique()
	samples_qc_pass = sample_qc.filter(qc_pass=True)['sample'].unique()
	samples = set(isescan_samples).intersection(set(samples_qc_pass))
	print(f'{len(isescan_samples)} ISESCan; {len(samples_qc_pass)} pass; {len(samples)} final samples')

	# Note: you could filter if est chrom cov is too low (like, below 10 or so) at this point
	contig_stats = estimate_relative_coverage(contig_stats)

	# Join per-contig ISEScan count info with per-contig assembly relative coverage info
	contig_ise = isescan_sum.join(contig_stats, on=['sample', 'contig']).filter(
		pl.col('family') != 'total', # This is implicitly done by the join-on contig, but nice to be explicit
		pl.col('sample').is_in(samples), # Filter to the QC-pass sample set
		pl.col('length') >= 700*pl.col('count') # Ignore contigs that are too small (too noisy / small piece of IS)
	)

	# Normalze IS couns by contig relative coverage
	contig_ise = contig_ise.with_columns(est_copy = pl.col('count')*pl.col('rel_cov'))

	# Define "common" IS families as ones where the average estimated copy per sample is at least 1
	common_fams = (
		contig_ise
		.group_by('family')
		.agg(pl.col('est_copy').sum() / len(samples)) # average count
		.filter(pl.col('est_copy') >= 1)
		.sort('est_copy', descending=True)
		.select('family').to_series().to_list()
	)
	# common_fams = contig_ise.group_by('family').agg(pl.col('count').sum()/len(samples)).sort('count', descending=True)
	# common_fams = ['IS1182', 'IS200/IS605', 'IS3', 'IS21', 'IS6', 'IS256', 'IS110']

	# Collapse uncommon fams into "other" category
	contig_ise = contig_ise.with_columns(
		family = pl.when(pl.col('family').is_in(common_fams)).then('family').otherwise(pl.lit('other'))
	)

	# Now, sum over sample (per family)
	sample_ise = contig_ise.group_by('sample', 'family').agg(
		count = pl.col('count').sum(),
		est_copy = pl.col('est_copy').sum()
	)

	# Make an "empty" dataframe with every combo of sample x family and then left-join to it,
	# and fill nulls with zeros -- so that we actually "count" zeros. 
	zero_ise = (
		pl.DataFrame(list(samples), schema=['sample'])
		.join(pl.DataFrame(common_fams+['other'], schema=['family']), how='cross')
	)
	sample_ise = zero_ise.join(sample_ise, how='left', on=['sample', 'family']).fill_null(strategy='zero')

	# Finally, join with sample metadata to get the year as well
	sample_ise = sample_ise.join(sample_qc['sample', 'year_disp'], on='sample')
	
	sample_ise.write_csv(output_csv)



def get_common_ise_families(isescan_sum):
	return (
		isescan_sum
		.filter(pl.col('family') != 'total')
		.filter(pl.col('dna_len') >= 700*pl.col('count'))
		.group_by('family')
		.agg(pl.col('count').sum())
		.sort('count', descending=True)
	)
	

# Heuristic to estimate the coverage of the chromosome, so that we can calculate
# the relative coverage of a given contig in the assembly. 
def estimate_relative_coverage(assembly_contig_stats):
	chromosome_coverage_estimate = (
		assembly_contig_stats
		.filter(pl.col('length') >= 5e3) # Only consider contigs >5k
		.sort(['sample', 'cov']) # Sort contigs by increasing coverage (per-sample) ...
		.with_columns(
			cum_len=pl.cum_sum('length').over('sample'), # ... so that you can get cumulative length
			total_length_g5k=pl.sum('length').over('sample') # + compute total assembly length (for long contigs)
		)
		.with_columns(frac = pl.col('cum_len') / pl.col('total_length_g5k'))
		
		# Finally, get the length-weighted median coverage by taking the coverage of
		# the first contig where the cumulative length is >0.5 the total length
		.filter(pl.col('frac') >= 0.5) # Note that `filter` preserves row order
		.unique(subset='sample', keep='first', maintain_order=True)
		.rename({'cov':'est_chrom_cov'})
		.select('sample', 'est_chrom_cov')
	)

	# Now, join this with the contig stats and compute relative coverage
	return (
		assembly_contig_stats
		.join(chromosome_coverage_estimate, on='sample')
		.with_columns(rel_cov = pl.col('cov')/pl.col('est_chrom_cov'))
	)




if __name__ == '__main__':
	main()