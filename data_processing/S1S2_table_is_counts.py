import polars as pl
import xlsxwriter # Not explicitly used, but need to verify this package is installed in order to write excel from pl.
import locate

# Collect all samples to be used in the collection.
ncbi = pl.read_csv('data/ncbi_all_filtered_samples_full.csv')
shc = pl.read_csv('data/shc_samples.csv').filter(pl.col('species') != 'undefined') # remove ERV_165
ente = (
	pl
	.read_csv('data/other_ente_samples.csv')
	# skip the reference genomes as they aren't really part of our collection
	.filter(collection=locate.SHC_ENTE)
	.with_columns(
		genus=pl.lit('Enterococcus'),
		# One is labeled 'durans-hirae'; it is durans according to GTDB-Tk.
		species=pl.col('sample').str.split('_').list.last().str.split('-').list.first()
	)
)
# Select only necessary columns
cols = ['sample', 'collection', 'genus', 'species']

df = pl.concat([ncbi.select(cols), shc.select(cols), ente.select(cols)])

# Simplify the collections to their base names to be more externally-friendly
coll = pl.col('collection')
df = df.with_columns(
	collection=pl
	.when(coll.str.starts_with(locate.NCBI_NANOPORE)).then(pl.lit(locate.NCBI_NANOPORE))
	.when(coll.str.starts_with(locate.NCBI_COMPLETE)).then(pl.lit(locate.NCBI_COMPLETE))
	.when(coll.str.starts_with(locate.SHC_ISO)).then(pl.lit(locate.SHC_ISO))
	.when(coll.str.starts_with(locate.SHC_ENTE)).then(pl.lit(locate.SHC_ENTE))
	.otherwise(locate.NA)
)

# Collection sizes
print(df['collection'].value_counts(sort=True))

# Count ISE per sample, for each IS family
ise = (
	pl
	.read_csv('data/isescan_per_contig.csv')
	# Count IS per sample per family
	.group_by('sample', 'family')
	.agg(pl.col('count').sum())
	# Pivot with sample as row idx, fam as cols
	.pivot(index='sample', on='family', values='count')
	.fill_null(0)
)

# Combine the df with the IS counts, sort rows/cols, and print out to file
df = (
	df
	# Note that SHC E. faecalis, and several E. faecium, did not complete ISEScan runs
	.join(ise[['sample'] + sorted(ise.columns[1:])], on='sample')
	.sort(by=['genus','species'])
)

# df.write_csv('data/supplementary_table_is_counts_per_sample.csv')

public = pl.col('collection').str.starts_with(locate.NCBI)
shc = pl.col('collection').str.starts_with(locate.SHC)

df.filter(public).write_excel('Table_S1_IS_counts_public_data.xlsx')
df.filter(shc).sort(by='collection').write_excel('Table_S2_IS_counts_SH_ente.xlsx')