import polars as pl
from ..config import data

# Here, we filter based on having geNomad & ISEScan results, and only species with >10 samples
df = pl.read_csv(data.ncbi.combined_qc_filtered)

# Remove "sp." species name, which has collapsed multiple species, except for Shigella, 
# which has intentionally not been further classified due to GTDB-tk ref
df = df.filter((pl.col('species') != 'sp.') | (pl.col('genus') == 'Shigella'))

# Import geNomad & ISEScan results to remove samples without those data
genomad = pl.scan_csv(data.genomad).select('sample').unique().collect().to_series().implode()
isescan = pl.scan_csv(data.isescan.summary).select('sample').unique().collect().to_series().implode()

df = (
	df
	.with_columns(
		genomad=pl.col('sample').is_in(genomad),
		isescan=pl.col('sample').is_in(isescan)
	)
	.filter(genomad=True, isescan=True)
	.drop('genomad', 'isescan')
)


# Only keep species with more than 10 samples at the very end
df = (
	df
	.join(df.group_by('genus','species').len(), on=['genus', 'species'], how='left')
	.filter(pl.col('len') > 10)
	.drop('len')
)

df.write_csv(data.ncbi.combined_final_samples)