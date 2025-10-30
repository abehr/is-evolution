import polars as pl
import xlsxwriter
from biobehr import io, gff_utils
from biobehr.gffpl import GFFPL
from efm.config import cfg, data

'''
Note that this requires additionally downloading annotation files (GFF3) from NCBI E. faecium.
The path to these can then be noted in the config.yaml.

Also, poppunk tree genomes list must be included as well. 
'''

# Samples used in PopPUNK tree
with open(cfg.source.poppunk_tree_genomes) as f:
	df = pl.DataFrame([l.strip() for l in f], schema=['sample'])

df = df.with_columns(
	gff=pl.format('{}/{}.gff3', pl.lit(str(cfg.source.gff)), 'sample')
)

df = io.filter_col_files_exist(df, 'gff')

gff_unmodified = GFFPL(df).valid_cds().lf.collect()


# Add info about proximal genes for each gene
gff = gff_unmodified.with_columns(
	prev=pl.col('product').shift(1).over('sample', 'seqid'),
	prev_strand=pl.col('strand').shift(1).over('sample', 'seqid'),
	prev_end=pl.col('end').shift(1).over('sample', 'seqid'),
	next=pl.col('product').shift(-1).over('sample', 'seqid'),
	next_strand=pl.col('strand').shift(-1).over('sample', 'seqid'),
	next_start=pl.col('start').shift(-1).over('sample', 'seqid')
)

# For new supplementary table:
# genes (except hypothetical protein and other transposase) with ISL3 (1) adjacent, and (2) upstream+cooriented
tpase = gff_utils.filter_IS_tpase(gff)
# Cap to the 15 most common (very rarely there is an incorrect name inference by the annotation software)
fams = tpase['fam'].value_counts(sort=True)[:15]['fam'].to_list()

is_sense = pl.col('strand') == '+'
next_dist = pl.col('next_start') - pl.col('end')
prev_dist = pl.col('start') - pl.col('prev_end')

tpase = (
	tpase
	# Add tpase names and tpase length
	.with_columns(
		tpase = gff_utils.extract_IS_name_from_annotation(fams),
		length = ((pl.col('end') - pl.col('start')+1)/3).cast(int)
	)
	# Establish upstream/downstream relative to the strandedness of the tpase, so we don't have to keep conditioning on it
	.with_columns(
		upstream_gene = pl.when(is_sense).then('prev').otherwise('next'),
		downstream_gene = pl.when(is_sense).then('next').otherwise('prev'),

		upstream_strand = pl.when(is_sense).then('prev_strand').otherwise('next_strand'),
		downstream_strand = pl.when(is_sense).then('next_strand').otherwise('prev_strand'),
		
		upstream_dist = pl.when(is_sense).then(prev_dist).otherwise(next_dist),
		downstream_dist = pl.when(is_sense).then(next_dist).otherwise(prev_dist)
	)
	.with_columns(
		upstream_coor = pl.col('strand') == pl.col('upstream_strand'),
		downstream_coor = pl.col('strand') == pl.col('downstream_strand'),
	)
	.select([
		'sample', 'seqid', 
		'upstream_gene', 'upstream_coor', 'upstream_dist',
		'fam', 'tpase', 'length', 
		'downstream_gene', 'downstream_coor', 'downstream_dist'
	])
)


# Specifically look into folT
folT = 'Folate transporter FolT'
folT_genomes = gff.filter(product=folT)['sample'].n_unique() # no. genomes with folT
folT_proximal = ((pl.col('upstream_gene') == folT) | (pl.col('downstream_gene') == folT))
isl3_folT_adjacent = tpase.filter(fam='ISL3').filter(folT_proximal)['sample'].n_unique()
isl3_folT_downstream_cooriented = tpase.filter(fam='ISL3', downstream_gene=folT, downstream_coor=True)['sample'].n_unique()
print(
	f"Of {folT_genomes} genomes with folT, {isl3_folT_adjacent} have an ISL3 directly adjacent to folT; " \
	f"{isl3_folT_downstream_cooriented} of these have ISL3 upstream of & co-oriented with folT."
)


# Generally, go through each IS family and see which genes are the most commonly adjacent, and the number of times they appear in each configuration. 
common_fams = ('ISL3', 'IS30', 'IS256', 'IS3', 'IS200/IS605', 'IS110', 'IS6')
df_fam = tpase.filter(pl.col('fam').is_in(common_fams))
upstream_gene_counts = (
	df_fam
	.group_by('fam', 'upstream_gene', 'upstream_coor').len().drop_nulls()
	.pivot(values='len', on='upstream_coor').fill_null(0)
	.rename(dict(upstream_gene='gene', true='upstream_coor', false='upstream_anti'))
)
downstream_gene_counts = (
	df_fam
	.group_by('fam', 'downstream_gene', 'downstream_coor').len().drop_nulls()
	.pivot(values='len', on='downstream_coor').fill_null(0)
	.rename(dict(downstream_gene='gene', true='downstream_coor', false='downstream_anti'))
)
gene_counts = upstream_gene_counts.join(downstream_gene_counts, on=['fam', 'gene'], how='full', coalesce=True).fill_null(0)
gene_counts = (
	gene_counts
	.with_columns(total=pl.sum_horizontal(gene_counts.select(pl.col(int))))
	.sort('total', descending=True)
)
common_genes = gene_counts.filter(
	~pl.col('gene').str.contains('(?i)transposase'), 
	pl.col('gene').str.to_lowercase() != 'hypothetical protein'
)

# For ISL3 only, also count the subset of IS→ X→ where the distance is in the potential promoter range,
# based on the common 3' flank length of ISL3. Tuned for this flank length; do not report for other 
# IS families, where this doesn't necessarily make sense.
downstream_coor_nearby_isl3 = (
	df_fam
	.filter(pl.col('downstream_coor'), pl.col('downstream_dist').is_between(75, 200))
	.group_by('fam', 'downstream_gene').len()
	.rename(dict(downstream_gene='gene', len='downstream_coor_nearby_isl3'))
	# .select('downstream_gene').to_series().alias('gene')
	# .value_counts(sort=True)
	# .rename(dict(count='downstream_coor_nearby_isl3'))
)
common_genes = common_genes.join(downstream_coor_nearby_isl3, how='left', on=['fam', 'gene']).fill_null(0)

# Naming schema is gene relative to IS, i.e. "downstream" means "gene downstream of IS"
output_columns = dict(
	fam = 'fam', # Don't rename this column; we will not use it in the output display.
	gene = 'Gene name',
	total = 'Total observations',
	downstream_coor = 'IS→ X→ (Co-oriented)',
	upstream_coor = 'X→ IS→ (Co-oriented)',
	downstream_anti = 'IS→ ←X (Convergent)', # gene is downstream of, and not co-oriented with, IS
	upstream_anti = '←X IS→ (Divergent)', # gene is upstream of, and not co-oriented with, IS
	downstream_coor_nearby_isl3 = 'IS→ [75-200bp] X→ (Co-oriented)' # Promoter-like activity possible from 3' flank of ISL3 given distance
)
output_col_order = list(output_columns.keys())

# Sort and rename columns for output
common_genes = common_genes.select(output_col_order).rename(output_columns)

with xlsxwriter.Workbook(data.output / 'Table_genes_with_adjacent_is.xlsx') as xl:
	bold = xl.add_format({"bold": True})
	for fam in common_fams: # Output in this order
		genes = common_genes.filter(fam=fam).drop('fam')

		if fam != 'ISL3': genes = genes.drop('IS→ [75-200bp] X→ (Co-oriented)')

		# Simpler output but then we can't format the column widths etc.:
		# genes[:500].write_excel(xl, xl.add_worksheet(fam.replace('/','_')))

		ws = xl.add_worksheet(fam.replace('/', '_'))
		
		for i,row in enumerate(genes.iter_rows()):
			if i >= 500: break # Only write the 500 most common to the output table
			for j,val in enumerate(row):
				ws.write(i+1, j, val) # start at row 1 so we don't overwrite header row
		
		for j,name in enumerate(genes.columns):
			ws.write(0, j, name, bold)
			ws.set_column(j, j, 50 if j == 0 else 18)

