import polars as pl
from pathlib import Path

join_on=['junction', 'sample']

test_data = Path('test_data/shc_st117')
pathstats = pl.read_json(test_data / 'pathstats.json')
joint_paths = pl.read_json(test_data / 'joint_paths.json')
junction_gff = pl.read_csv(test_data / 'junction_gff.tsv')

# Analyze ISL3
jcts_to_study = pathstats.filter(tpase_only=True, isfam='ISL3') #, npaths=2)

# Analyze other common IS families
jcts_to_study = pathstats.filter(
	pl.col('tpase_only'),
	pl.col('isfam').is_in(['IS30','IS256','IS110','IS3','IS6']),
)


dg_gff = (
	# Filter joint_paths to the set of simple binary junctions we're looking at
	joint_paths.filter(pl.col('junction').is_in(jcts_to_study['junction']))
	# Join with the junction GFF and add a few helper columns
	.join(junction_gff.rename(dict(seqid='junction')), on=join_on, how='left', maintain_order='left_right')
	.with_columns(gene_length=pl.col('end')-pl.col('start'))
	.with_columns(gene=pl.struct('gene_length','product'))
	.with_columns(gene_pos = pl.struct('start','end', 'strand', 'product'))
	.with_columns(acc_l=pl.col('position').list.get(1), acc_r=pl.col('position').list.get(-2))
)

# Filter these junction gffs to the accessory region PLUS the genes that are immediately before/after, 
# but *completely* absolutely outside of, the accessory region. Basically we want to look at the accessory region
# plus slightly outside of it, to make sure that the backbone is correct and take a wider view,
# but we don't want to look across the entire flanking core blocks which can contain SNPs/indels
end_before = pl.col('end') < pl.col('acc_l')
start_after = pl.col('start') > pl.col('acc_r')
gene_within_acc = (pl.col('start') >= pl.col('acc_l')) & (pl.col('end') <= pl.col('acc_r'))

# True if the gene is completely or partially within accessory, or whether accessory is within gene.
# Also true if the gene touches the accessory right at the boundary. 
gene_overlaps_acc = (pl.col('end') >= pl.col('acc_l')) & (pl.col('start') <= pl.col('acc_r'))

dg_gff = (
	dg_gff
	.join(dg_gff.filter(end_before).group_by(join_on).agg(l0=pl.col('end').max()), how='left', on=join_on)
	.join(dg_gff.filter(start_after).group_by(join_on).agg(r0=pl.col('start').min()), how='left', on=join_on)
	.with_columns(l0=pl.col('l0').fill_null(strategy='one'), r0=pl.col('r0').fill_null(strategy='max'))
	.filter(pl.col('end') >= pl.col('l0'), pl.col('start') <= pl.col('r0'))
)

# Now, collapse the gffs to list of genes. We do this in two different ways for each junction:
#  - "accessory-cured" (ordered list of genes, skipping any that are completely within accessory)
#  - full list of genes as well.
path_grp = ['junction', 'sample', 'path_str']
dgg_acc_cured = (
	dg_gff.filter(~gene_within_acc)
	.group_by(path_grp, maintain_order=True)
	.agg(genes_cured = pl.col('gene'))
)
dgg_overlap_acc = (
	dg_gff.filter(gene_overlaps_acc & ~gene_within_acc)
	.group_by(path_grp, maintain_order=True)
	.agg(genes_overlap_acc = pl.col('gene'))
)
dgg = (
	dg_gff
	.group_by(path_grp, maintain_order=True)
	.agg(genes = pl.col('gene'), genes_pos = pl.col('gene_pos')) # all genes in the region, incl within accessory
	.join(dgg_acc_cured, on=path_grp, how='left')
	.join(dgg_overlap_acc, on=path_grp, how='left')
)

# An empty path has a possible insertional mutagenesis if it has a gene 
# that overlaps with where the accessory region is inserted (in the filled path).
dis = (
	dgg
	.with_columns(path_len = pl.col('path_str').str.count_matches('_')+1)
	.with_columns(
		possibly_insertionally_mutagenized=pl
		.when(path_len=2)
		.then(pl.col('genes_overlap_acc').fill_null([]).list.len() > 0)
		.otherwise(pl.lit(False))
	)
)

genes_disrupted = []
for j in dis.filter(possibly_insertionally_mutagenized=True)['junction'].unique():
	core_genes_overlapped = dis.filter(junction=j, possibly_insertionally_mutagenized=True)['genes_overlap_acc'].explode().unique()
	acc_genes = dis.filter(pl.col('junction') == j, pl.col('path_len') > 2)['genes'].explode().unique()
	missing = [x for x in core_genes_overlapped if x not in acc_genes]
	genes_disrupted.append(dict(
		junction = j,
		core_genes_overlapped = core_genes_overlapped.to_list(),
		unique_acc_genes = acc_genes.to_list(),
		disrupted = missing,
		is_disrupted = len(missing) > 0
	))
genes_disrupted = pl.DataFrame(genes_disrupted)


print(f'---- Genes disrupted: {len(genes_disrupted.filter(is_disrupted=True))} / {len(jcts_to_study)} junctions')
print(genes_disrupted.filter(is_disrupted=True).explode('disrupted').unnest('disrupted')['product'].sort().to_list())


###################################
# Distance away from nearest gene
# should only apply to intergenic ones
# NOTE: when a junction has 2 genes in its accessory region, it can be difficult to determine which one is the "new" insertion. 
# In that case we'll just pick the furthest one. And if they're in opposite orientations we'll call it ambiguous,
# and we will not consider it, rather than assuming incorrectly. 
###################################

accessory_ise = (
	dg_gff
	.filter(pl.col('path_length') > 2, gene_within_acc)
	.with_columns(acc=pl.struct('acc_l', 'acc_r', 'start', 'end', 'strand', 'product'))
	.select('junction', 'sample', 'path_str', 'path_length', 'acc')
)

ambiguous = []
distance = []
for (junction,sample),acc in accessory_ise.group_by(join_on):
	# ignore cases where we know a gene disruption occurred. 
	if junction in genes_disrupted.filter(is_disrupted=True)['junction']: continue
	# If there are multiple genes in the accessory, pick the more-downstream one.
	# If they have opposite orientations, it's too ambiguous and we can't decide anything about it. 
	strand = acc.unnest('acc')['strand'].unique()
	if len(strand) > 1:
		ambiguous.append(junction)
		continue

	r = acc.row((-1 if strand.item() == '+' else 0), named=True)

	g = (
		dg_gff
		.filter(junction=r['junction'], sample=r['sample'], path_str=r['path_str'])
		.select('start', 'end', 'strand', 'product')
	)
	if r['acc']['strand'] == '+':
		next_gene = g.filter(pl.col('start') >= r['acc']['end']).row(0, named=True)
		dist_from_ise = next_gene['start'] - r['acc']['end']
		dist_from_acc = next_gene['start'] - r['acc']['acc_r']
	else:
		next_gene = g.filter(pl.col('end') <= r['acc']['start']).row(-1, named=True)
		dist_from_ise = r['acc']['start'] - next_gene['end']
		dist_from_acc = r['acc']['acc_l'] - next_gene['end']
	
	# r['downstream'] = next_gene
	r['ise'] = r['acc']['product']
	r['d_dist_from_ise'] = dist_from_ise
	r['d_dist_from_acc'] = dist_from_acc
	r['d_strand_parity'] = r['acc']['strand'] == next_gene['strand']
	r['d_product'] = next_gene['product']
	distance.append(r)

distance = pl.DataFrame(distance).drop('acc')

x = (
	distance
	.group_by(['junction', 'path_str'])
	.agg(
		ise = pl.col('ise').unique(),
		gene = pl.col('d_product').unique(),
		ise_dist = pl.col('d_dist_from_ise').mean(),
		ise_std = pl.col('d_dist_from_ise').std(),
		acc_dist = pl.col('d_dist_from_acc').mean(),
		acc_std = pl.col('d_dist_from_acc').std(),
		d_strand_parity = pl
			.when(pl.col('d_strand_parity').n_unique() == 1)
			.then(pl.col('d_strand_parity').first())
			.otherwise(pl.lit(None)),
	)
	.fill_null(0, matches_supertype=True) # fill null standard deviations
	.filter(pl.col('d_strand_parity').is_not_null()) # shouldn't be any of these anyway
	.with_columns(
		ise=pl.col('ise').list.join(', '),
		gene=pl.col('gene').list.join(', ') # rare, if downstream gene varies across samples even within the same path in the same junction.
	) 
)

# Note that this does not exactly equal the number of junctions.
# We filter out ambiguous cases, but also if a junction has multiple possible "filled"
# paths, i.e. not binary, these are considered distinct. 
x['ise_dist'].mean() # ISL3: 183 nt (sd 110) // other: 286 nt (sd 191)
x['d_strand_parity'].value_counts(sort=True) # ISL3: 93 opposite, 8 parallel // Other: 18 opposite, 13 parallel

x.write_ndjson('ISE_intergenic_distance.jsonl')