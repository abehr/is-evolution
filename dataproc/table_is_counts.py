import polars as pl
import xlsxwriter # Not explicitly used, but need to verify this package is installed in order to write excel from pl.
from efm.config import data
from efm import dataset_definitions as defs

# Collect all samples to be used in the collection.
ncbi = pl.read_csv(data.ncbi.combined_final_samples)
shc = pl.read_csv(data.shc.iso)
ente = (
	pl
	.read_csv(data.shc.other_ente)
	# skip the reference genomes as they aren't really part of our collection
	.filter(collection=defs.SHC_ENTE)
	.with_columns(
		genus=pl.lit('Enterococcus'),
		# One is labeled 'durans-hirae by genome name'; it is durans according to GTDB-Tk.
		species=pl.col('sample').str.split('_').list.last().str.split('-').list.first()
	)
)

# Grab # contigs & assembly length from the assembly info file for the atypical ente (which is not already baked in)
ente_assembly_info = (
	pl.read_csv(data.assembly_info)
	.filter(pl.col('sample').is_in(ente['sample'].implode()))
	.select('sample', 'contigs', 'genome_length')
)
ente = ente.join(ente_assembly_info, on='sample')

# Select only necessary columns
cols = ['sample', 'collection', 'genus', 'species', 'genome_length', 'contigs']

df = pl.concat([ncbi.select(cols), shc.select(cols), ente.select(cols)])

# Simplify the collection names to be more externally friendly and match the manuscript references
df = df.with_columns(collection=defs.convert_collection_name)

# Collection sizes
print(df['collection'].value_counts(sort=True))

# Count ISE per sample, for each IS family
ise = (
	pl
	.read_csv(data.isescan.per_contig)
	# Count IS per sample per family
	.group_by('sample', 'family')
	.agg(pl.col('count').sum())
	# Pivot with sample as row idx, fam as cols
	.pivot(index='sample', on='family', values='count')
	.fill_null(0)
)

# Combine the df with the IS counts, sort rows/cols, and print out to file
ise_col_order = sorted(ise.columns[1:])
df = (
	df
	# Note that SHC E. faecalis, and several E. faecium, did not complete ISEScan runs
	.join(ise[['sample'] + ise_col_order], on='sample')
	.sort(by=['genus','species'])
)

# df.write_csv('data/supplementary_table_is_counts_per_sample.csv')

public = pl.col('collection').str.starts_with(defs.C_NCBI)
shc = pl.col('collection').str.starts_with(defs.C_SHC)

# Finally, get contig categories per sample with geNomad for the public data
# (we didn't analyze SH data with geNomad)
genomad_col_order = ['chromosome', 'plasmid', 'virus', 'undetermined']
genomad_contig_types = (
	pl
	.read_csv(data.genomad)
	.group_by('sample', 'contig_type').len()
	.pivot(index='sample', on='contig_type').fill_null(0)
)


# Write out public data to Table S1, including geNomad contig classifications
(
	df
	.filter(public)
	.join(genomad_contig_types, on='sample')
	.select(cols + genomad_col_order + ise_col_order)
	.write_excel(data.output / 'Table_IS_counts_public_data.xlsx')
)

# Write out SHC data to Table S2 (no geNomad info)
df.filter(shc).sort(by='collection').write_excel(data.output / 'Table_IS_counts_SH_ente.xlsx')