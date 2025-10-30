# Nextflow pipelines to support timeline of historical IS counts from short read data

This Nextflow pipeline is split into two workflows. Example configuration & params files are provided.

1. `pre.nf` takes in a set of SRA run accessions (Illumina paired-end isolate data), which is produced by the NCBI entrez fetch/preprocessing/filter steps. Reads are downloaded, deduplicated, trimmed, sketched, QC-filtered, downsampled (if necessary), and assembled (if estimated coverage is sufficient).
2. `post.nf` takes in these assemblies and runs ISEScan, FastANI, and CheckM2.

Additional executor and process configuration, including pipeline dependencies (e.g. by specifying conda envs or Docker containers), must be defined in a `nextflow.config` in order to run the workflow. The dependencies are shown in `env_nf_sr.yml`. However, we recommend configuring a one-off conda environment separately for the CheckM2 process, so that you can initialize it and configure the databases. 

After completing the pipeline, initial post-processing to concatenate the results from the different pipeline steps can look something like this:

```bash
cd nf_result

# FastANI
cat fastANI/*.txt > fastANI.tsv

# Mash
cat mash_dist/*.txt > mash_dist_reads.tsv

# CheckM2 (strip header row except one)
awk 'FNR==1 && NR!=1 { next } { print }' checkM2/*.txt > checkM2.tsv

# SPAdes: get contig headers paired with sample ID
for f in assembly/*.fna; do
	stem=$(basename "$f" .fna)
	grep '^>' "$f" | awk -v fname="$stem" '{print fname "\t" $0}'
done > assembly_contig_headers.tsv
```

```python
import polars as pl
from pathlib import Path
from ...data_processing.isescan import fast_parse_isescan_sum_files
from ...common import io

# SPAdes: parse info from contig headers
spa = (
	pl
	.read_csv('assembly_contig_headers.tsv', separator='\t', has_header=False, new_columns=['sample', 'contig'])
	.with_columns(
		sample=pl.col('sample').str.strip_suffix('.fna'),
		contig=pl.col('contig').str.strip_prefix('>')
	)
	.with_columns(data=pl.col('contig').str.split_exact('_', n=5))
	.with_columns(
		length=pl.col('data').struct[3].cast(int),
		cov=pl.col('data').struct[5].cast(float)
	)
	.drop('data')
)
spa.write_csv('assembly_contigs_parsed.csv')


# ISEScan
ise_sum_files = io.fast_glob(Path('./isescan'), '*.sum')
ise = (
	fast_parse_isescan_sum_files(ise_sum_files)
	# .with_columns(pl.col('sample').str.strip_suffix('.fna')) # Not necessary actually for the efm one ?? 
	# .filter(pl.col('family') != 'total') # get rid of sum total count line
)
ise.write_csv('isescan_sum.csv')
```

