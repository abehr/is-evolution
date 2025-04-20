import polars as pl
import locate

c = pl.read_csv('data/ncbi_complete_datasets.csv')
n = pl.read_csv('data/ncbi_nanopore_samples.csv')

# Parse organism column; ignore subspecies here for parity with GTDB-Tk result
c = c.with_columns(
	genus=pl.col('organism').str.split(' ').list.get(0),
	# This method fills null if the species name has more than 2 words in it.
	# species=pl.col('organism').str.extract(r'^\S+\s+(\S+)$').fill_null('')
	species=pl.col('organism').str.split(' ').list.get(1).fill_null('')
)

# Parse GTDB-tk taxonomy; ignore sub-genus here for parity with ncbi complete
n = n.with_columns(
	# Ignore the sub-genus designation that GTDB-Tk can give (e.g. "Enterococcus_B") for parity with NCBI txid type
	genus=pl.col('genus').str.split('_').list.get(0),
	# Here, we don't have to worry about whether the species name has >2 words
	# so we can use a regex.
	# but sometimes it's genus-only, which means that splitting it into a list
	# and taking the 2nd item would cause an error. 
	species=pl.col('species').str.extract(r'^\S+\s+(\S+)$').fill_null('')
)

# Finally, we also need to ignore the subspecies clade (e.g. "hormaechei_B") for parity with NCBI txid type
n = n.with_columns(species=pl.col('species').str.split('_').list.get(0))

# Compute average genome length per genus
c_genome_lengths = (
	c.filter(pl.col('genome_fp').is_not_null(), ~pl.col('dupe'), pl.col('status') == 'current')
	.select('genus', 'genome_length')
)
n_genome_lengths = (
	n.filter(pl.col('genome_fp').is_not_null(), pl.col('genus').is_not_null())
	.select('genus', 'genome_length')
)
mean_length = (
	pl
	.concat([c_genome_lengths, n_genome_lengths])
	.group_by('genus')
	.agg(
		nsamp=pl.len(),
		mean_length=pl.mean('genome_length'),
		std_length=pl.std('genome_length')
	)
	.filter(pl.col('nsamp') >= 10)
	.sort(by='genus')
)

mean_length = mean_length.with_columns(
	min_cutoff = pl.col('mean_length') - 2*pl.col('std_length'),
	max_cutoff = pl.col('mean_length') + 2*pl.col('std_length'),
)

# ==========================
# Filter to good assemblies. This means slightly different things for the complete vs nanopore.
# ==========================

c_filtered = (
	c.join(mean_length, how='left', on='genus')
	.filter(
		~pl.col('dupe'),
		# pl.col('status') == 'current',
		pl.col('genome_fp').is_not_null(),
		pl.col('genus').is_in(locate.genera),
		pl.col('species') != '', # not empty species either
		pl.col('genome_length').is_between('min_cutoff', 'max_cutoff'),
		pl.col('contigs') <= 20
	)
)

n_filtered = (
	n.join(mean_length, how='left', on='genus')
	.filter(
		pl.col('genome_fp').is_not_null(),
		pl.col('genus').is_in(locate.genera),
		pl.col('species').is_not_null(),
		pl.col('species') != '', # not empty species either
		pl.col('genome_length').is_between('min_cutoff', 'max_cutoff'),
		pl.col('contigs') <= 20,
		pl.col('mean_coverage') >= 40
	)
)

# Combine all into a single dataset
cols = ['sample', 'collection', 'genus', 'species', 'assembly', 'genome_length', 'contigs']
all_filtered = pl.concat([
	c_filtered.select(cols),
	n_filtered.with_columns(assembly=pl.lit(None)).select(cols)
])

assert all_filtered['sample'].n_unique() == len(all_filtered), 'Sample ID collision'

all_filtered = all_filtered.with_columns(eskape=locate.eskape)

all_filtered.write_csv('data/ncbi_all_filtered_samples.csv')