import polars as pl
from pathlib import Path
from ..common.dataset_definitions import survey_genera, NCBI_COMPLETE
from ..common.io import col_exists
from ..config import cfg, data


cols = dict(
	accession='sample',
	pairedAccession='paired_accession',
	sourceDatabase='source',
	# -- assembly info:
	assemblyName='assembly',
	assemblyStatus='status',
	# note you could further unnest biosample and get sample accession & collection date
	# collectionDate='collection_date',
	# -- assembly stats:
	numberOfContigs='contigs',
	totalSequenceLength='genome_length',
	# -- organism:
	organismName='organism'
	# note you could unnest ANI and use additional info to determine species
)


# We have a taxon dataset for each genus
def get_ncbi_complete(base_fp, taxon_name):
	fp = base_fp / taxon_name / 'ncbi_dataset/data'
	return (
		pl
		.scan_ndjson(fp / 'assembly_data_report.jsonl')
		.unnest('assemblyInfo', 'assemblyStats', 'organism')
		.select(list(cols.keys()))
		.rename(cols)
		.with_columns(collection=pl.lit(f'{NCBI_COMPLETE}_{taxon_name}'))
		# synonymous with locate.ncbi_complete_genome & complete_genome_name
		.with_columns(genome_fp=pl.format(
			'{}/{}/{}_{}_genomic.fna',
			pl.lit(str(fp)), 'sample', 'sample', 
			pl.col('assembly').str.replace_all(' ', '_')
		))
	)


def polish_ncbi_complete(df):
	print('Resolve genome files exist')
	df = df.with_columns(col_exists('genome_fp'))

	print('Resolve duplicate accessions')
	unique_accessions = (
		df
		.with_columns(
			prefer_exists = pl.when(pl.col('genome_fp').is_not_null()).then(0).otherwise(1),
			prefer_refseq = pl.when(source='SOURCE_DATABASE_REFSEQ').then(0).otherwise(1)
			# prefer_current = pl.when(status='current').then(0).otherwise(1),
		)
		.sort(by=['prefer_exists', 'prefer_refseq']) # 'prefer_current', 'prefer_refseq'])
		.unique(subset='assembly', keep='first')
		.select('sample')
		.to_series()
		.implode()
	)

	# Note which accession should be kept, but do not remove the duplicates from the df
	df = df.with_columns(dupe=~pl.col('sample').is_in(unique_accessions))

	return df



# Ultimately, should only use assemblies that are dupe=False, assembly exists, and status='current'
if __name__ == '__main__':
	base_fp = Path(cfg.source.ncbi_complete_datasets)
	print(f'Load NCBI complete metadata for {len(survey_genera)} genera')
	lf = pl.concat([get_ncbi_complete(base_fp, g.lower()) for g in survey_genera])
	# Note that genome_length col comes in as string, but we don't need to cast it for now. 
	# df = df.with_columns(pl.col('genome_length').cast(int))
	print('Collect DataFrame')
	df = lf.collect()
	
	df = polish_ncbi_complete(df) 
	df.write_csv(data.ncbi.complete_genomes)
