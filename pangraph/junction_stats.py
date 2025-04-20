from pathlib import Path
import polars as pl
import sys
sys.path.append('..')
from util import io
from util.gffpl import GFFPL

from subgraph_utils import parse_joints, joints_pos_with_spacers, extract_junction_gff, parse_subgraph_paths
from junction_stat_utils import *

test_data = Path('test_data/shc_117_efm')
save_objects = False
if len(sys.argv) > 1:
	test_data = Path(sys.argv[1])
	assert test_data.is_dir(), f'Path {test_data} is not a dir'

# 1. get the junction positions & convert the GFFs to get them per-sample-junction
print('1. Load GFFs & junction positional info & genome fasta files')
joints_pos = parse_joints(test_data / 'backbone_joints/joints_pos.json')

# For every dataset, we should verify here that all three of the spacer column metrics agree,
# otherwise the spacers are not neatly confined to accessory regions and there's something unexpected.
joints_pos = joints_pos_with_spacers(joints_pos, test_data / 'input_fasta')
assert len({joints_pos[x].sum() for x in ('jct_contains_N', 'n_spacers_jct', 'n_spacers_acc')}) == 1

# Dump junction gff to file so we don't have to reload it every time, as it is the slowest step.
print('2. Load junction GFFs')
junction_gff_fp = test_data / 'junction_gff.tsv'
if not junction_gff_fp.is_file():
	gff_files = io.fast_glob(test_data / 'annotations_gffpl_preformatted', '*', False)
	gff = {x.stem:GFFPL(x).cds_rna().lf.collect() for x in gff_files}
	extract_junction_gff(gff, joints_pos).write_csv(junction_gff_fp, separator='\t')
else:
	print('--- Using precomputed junction GFFs')

junction_gff = pl.read_csv(junction_gff_fp, separator='\t')

# 2. Get the subgraphs
print('\n3. Load subgraph paths; generate per-junction path statistics')
joint_paths = parse_subgraph_paths(test_data / 'backbone_joints/joints_pangraph')

# If you want to filter to a subset of the samples, you can do that here. 
# joint_paths = joint_paths.filter(pl.col('sample').str.starts_with('EF') & (pl.col('sample') != 'EF_B_17'))

# Remove artifacts from stitched contigs.
# Note that filtering to n_spacers_acc=0 also filters out nulls, but that's ok to merge with joint_paths
joint_paths = joint_paths.join(joints_pos.filter(n_spacers_acc=0)['junction','sample'],  on=['junction','sample'])

# Let's also make sure that our subgraphs look right.
untrue_junctions_disagreement = (
	joint_paths
	.filter(pl.col('path_length') == 1) # junctions that aren't true junctions
	.group_by('junction')
	.agg(path_str=pl.col('path_str').n_unique())
	.filter(pl.col('path_str') != 1) # there should only be one possible path through these
)
empty_junctions_disagreement = (
	joint_paths
	.filter(pl.col('path_length') == 2) # empty junctions
	.group_by('junction')
	.agg(path_str=pl.col('path_str').n_unique())
	.filter(pl.col('path_str') != 1) # there should only be one possible path through these
)
filled_junctions_disagreement = (
	joint_paths
	.filter(pl.col('path_length') >= 2) # filled junctions
	.with_columns(
		first_block = pl.col('path').list.first(),
		last_block = pl.col('path').list.last()
	)
	.group_by('junction')
	.agg(
		first_block = pl.col('first_block').n_unique(),
		last_block = pl.col('last_block').n_unique()
	)
	# the first and last block of the junction should usually correspond to the flanking core blocks.
	.filter((pl.col('first_block') != 1) | (pl.col('last_block') != 1))
)
assert len(untrue_junctions_disagreement) == 0
assert len(empty_junctions_disagreement) == 0
print(f'Found {len(filled_junctions_disagreement)} junctions where edge core blocks did not always resolve the same')
# if this is >0, it's likely fine, likely due to repetitive sequence between core & accessory. just worth noting. 
# it's just worth noting, because we are less able to distinguish what is truly the "accessory" region of these blocks, in these cases.

pathstats = junction_path_stats(joint_paths)

print('\n4. Define junction genotypes')
accessory_gff = junction_accessory_gff(joint_paths, junction_gff)

# TODO: redo a bit more rigorously; make sure understanding of this file type is correct. 
exclude_subset = pl.col('in_acc').list.set_difference(['subset'])
genomad_viral = io.pl_read(test_data / 'annotations/junct_pos/genomad_real.csv').collect()
genomad_viral = (
	genomad_viral
	.rename({'iso':'sample'})
	.filter(in_core='subset')
	.with_columns(viral_length=pl.col('ie')-pl.col('ib')) # may need to account for wrapping if we use this
	.group_by(['junction', 'sample'])
	# Genomad often finds multiple viral regions within a single junction, but aggregate so we have only one.
	.agg(
		num_viral = pl.len(),
		viral_length = pl.col('viral_length').sum(),
		in_acc = pl.col('in_acc').unique()
	)
	# If there are multiple viral regions, and one is 'start' and one is 'subset',
	# the consensus is basically 'start' (bc overall, the viral regions are not fully contained)
	# in the accessory and they overlap start). If we find 'start' and 'end', it should be 'superset',
	# although I don't find any of those. Sort of basing off of https://github.com/mmolari/ecoliST131-structural-evo/blob/2d41d70a155c8753a42e2f28e908dd4654c8d68c/scripts/annotations/assing_junction.py#L9
	.with_columns(
		in_acc=pl.when(exclude_subset.list.len() == 0).then(pl.lit('subset'))
		  .when(exclude_subset.list.len() == 1).then(exclude_subset.list.first())
		  .otherwise(pl.lit('superset'))
	)
)

'''
Here we can quickly see ALL of the genes in accessory across all isolates. 
This is also really helpful bc we can look at what's hiding in "other" 
and decide if we want to make another category. 
Note that IS3 count would be slightly inflated bc of orfA/orfB.

NOTE: we could also filter on 'family', not just 'transposase',
because this seems to filter to IS elements more strictly. 
'''
accessory_gff.select('isfam').to_series().value_counts(sort=True)

accessory_paths = (
	accessory_gff

	# Group & aggregate to get the list of gene products within each accessory region
	.group_by(['junction', 'sample', 'path_str', 'path_length'], maintain_order=True)
	.agg(products='product', isfams='isfam')

	# If # tpase = # genes, then this is a tpase-only path
	.with_columns(
		gene_count = pl.col('products').list.len(),
		tpase_count = pl.col('isfams').list.drop_nulls().list.len(),
		isl3_count = pl.col('isfams').list.eval(pl.element() == 'ISL3').list.sum()
	)
	.with_columns(tpase_only=(pl.col('tpase_count') == pl.col('gene_count')))

	# If the first or last gene has a named tpase fam, then it is an IS-bounded path
	.with_columns(
		ise_bounded=(
			(pl.col('isfams').list.first().is_not_null() & (pl.col('isfams').list.first() != 'none')) |
			(pl.col('isfams').list.last().is_not_null() & (pl.col('isfams').list.last() != 'none'))
		)
	)

	# TODO: also consider the gene at the edge of the core region...

	# For IS-bounded paths, only look at first & last annotation for consensus IS fam
	.with_columns(
		isfams=pl
		.when(ise_bounded=True)
		.then(pl.col('isfams').list.gather([0,-1]))
		.otherwise('isfams')
	)
	
	# Decide on a consensus IS family annotation (if there is one).
	# Note that isfam is kinda not valid if tpase_only and ise_bounded are both false. 
	# In that case, it just tells you what the consensus IS familiy is, if the accessory path has them. 
	.with_columns(
		isfam=pl.col('isfams')
				.list.unique()
				.list.set_difference(['none', None])
	)
	.with_columns(
		isfam=pl
		.when(pl.col('isfam').list.len() == 1).then(pl.col('isfam').list.first())
		.when(pl.col('isfam').list.len() > 1).then(pl.lit('multi'))
		.otherwise(pl.lit(None))
	)
	.drop('isfams')

	# Add genomad viral region info
	.join(genomad_viral, on=['junction', 'sample'], how='left')
	.fill_null(0) # note that 'in_acc' column is still
	.with_columns(
		accessory_phage=pl
		# Optionally, a more nuanced view based on partial overlap with accessory region.
		
		# .when((pl.col('num_viral') > 0 ) & (pl.col('in_acc') == 'subset')).then(pl.lit('yes'))
		# .when((pl.col('num_viral') > 0 ) & (pl.col('in_acc') != 'subset')).then(pl.lit('partial'))
		# .when(num_viral=0).then(pl.lit('no'))
		# .otherwise(pl.lit('unexpected'))

		# I think the better thing to do w.r.t. this is reasses within the subgraph by updating the coordinates. For now, simplifying. 
		.when(num_viral=0).then(pl.lit(False)).otherwise(pl.lit(True))
	)
)

# This utilizes the trinary-ness of polars boolean logic. e.g. if some paths within the junction are 
# tpase-only, but others are not, then we cannot definitively say that the junction is 
# tpase-only or not tpase-only. Thus we assign it a null value for tpase-only. "ambiguous" would be more accurate. 
resolve_ambiguous = lambda c: pl.when(pl.col(c).list.len() == 1).then(pl.col(c).list.first()).otherwise(pl.lit(None))


accessory_pathstats = (
	accessory_paths
	.group_by('junction')
	.agg(
		n_accessory_paths=pl.col('path_str').n_unique(), # should correspond with number of non-empty paths
		n_accessory_samp=pl.col('sample').len(),
		accessory_pangenome_length=pl.col('gene_count').sum(), # this is not nonredundant, so it is kind of a misnomer
		n_accessory_tpase=pl.col('tpase_count').sum(),
		n_accessory_isl3=pl.col('isl3_count').sum(),
		tpase_only=pl.col('tpase_only').unique(),
		ise_bounded=pl.col('ise_bounded').unique(),
		isfam=pl.col('isfam').unique(),
		n_accessory_samp_overlapping_phage=pl.col('accessory_phage').sum()
	)
	.with_columns(
		tpase_only=resolve_ambiguous('tpase_only'),
		ise_bounded=resolve_ambiguous('ise_bounded')	
	)
	# Should abstract this so that we are not copy-pasting this twice.
	.with_columns(isfam=pl.col('isfam').list.unique().list.set_difference(['none',None]))
	.with_columns(
		isfam=pl
		.when(pl.col('isfam').list.len() == 1).then(pl.col('isfam').list.first())
		.when(pl.col('isfam').list.len() > 1).then(pl.lit('multi'))
		.otherwise(pl.lit(None))
	)

	# a given junction is definitively phage-like if every accessory path overlaps with a phage;
	# if some of the paths do then it's ambiguous; otherwise it's non-phage. 
	.with_columns(
		is_phage=pl
		.when(pl.col('n_accessory_samp_overlapping_phage') == pl.col('n_accessory_samp')).then(pl.lit('yes'))
		.when(pl.col('n_accessory_samp_overlapping_phage') == 0).then(pl.lit('no'))
		.otherwise(pl.lit('ambiguous'))
	)
)

'''
Note that now, things can be "null" for tpase_only (e.g.) even if they were not
evaluated that way above -- all of the accessory_pathstats cols will be null 
for any junction where there is only 1 path (not a true junction), or where
no path has any accessory genome (e.g. small accessory path with no annotation).
I would perhaps consider those also sort of not true paths. 
'''
pathstats = pathstats.join(accessory_pathstats, on='junction', how='left')

'''
# Side-quest, I guess: it's probably worth looking at some of these ambiguous junctions
ambiguous_ise_junctions = pathstats.filter(is_junction & pl.col('n_accessory_paths').is_not_null() & pl.col('tpase_only').is_null())
# For these junctions, which paths are not ISE-only? What do they look like? 
interesting_paths = accessory_paths.filter(pl.col('junction').is_in(ambiguous_ise_junctions) & ~pl.col('tpase_only'))
for x in interesting_paths.select('products').unique().to_series().to_list():print('\n'.join(x) + '\n')
# TODO: these are very interesting and have tons of IS elements in them. These may actually be the most hotspot-prone areas, given that they have diverse paths. 
'''

# be careful with boolean logic because of the trinary thing with bool type that pl does.
true_junction = pl.col('npaths') > 1
has_accessory = pl.col('n_accessory_paths') > 0
ise_only = (pl.col('tpase_only') == True) # at this point it's only ISE because of the family restriction (no 'none')
ise_bounded = (pl.col('tpase_only').is_null() | (pl.col('tpase_only') == False)) & (pl.col('ise_bounded') == True)

j = pathstats.filter(true_junction)
a = j.filter(has_accessory)
ise = a.filter(ise_only)
other = a.filter(~pl.col('junction').is_in(ise.select('junction')))
phage = other.filter(pl.col('is_phage').is_in(['yes','ambiguous']))


skd = {
	'root': {'Junctions': len(j)},
	'Junctions': {
		# 'Empty junctions': len(j)-len(a),
		'Gene-containing junctions': len(a)
	},
	'Gene-containing junctions': {
		'IS-only activity': len(ise),
		# 'Complex activity': len(a)-len(ise)
	},
	# 'Complex activity': {
		# 'Phage-associated': len(phage),
		# 'Other': len(other)-len(phage)
	# },
	'IS-only activity': {r[0]:r[1] for r in ise.select('isfam').to_series().value_counts().iter_rows()}
}


fig_file = Path(test_data / 'sankey.pdf')
print(f'Writing image to file at {fig_file.absolute()}')
with open(test_data / 'sankey.txt', 'w') as f: f.write(str(skd))
fig = sankey(skd)
fig.write_image(fig_file)
fig.show()



############################################
# Helper functions just for interactive work
############################################

def view_junction_gff(junction, sample):
	path_gff = junction_gff.filter(seqid=junction, sample=sample)
	joint = joint_paths.filter(junction=junction, sample=sample).row(0, named=True)

	for i,block in enumerate(joint['path']):
		block_start = joint['position'][i]
		block_end = joint['position'][i+1]
		print(f'--- {block} ({block_start}-{block_end})')
		for g in path_gff.filter(pl.col('start').is_between(block_start, block_end, closed='left')).iter_rows(named=True):
			print(f"{g['start']}-{g['end']} ({g['strand']})\t{g['product']}")

def view_pathstats(df):
	cols = ['nsamp', 'n_accessory_samp', 'npaths', 'n_accessory_paths', 'accessory_pangenome_length', 'n_accessory_tpase', 'n_accessory_isl3']
	print("JUNCTION\t\t\tsamples\tsamples\tpaths\tpaths\tgenes\ttpase\tISL3")
	for r in df.iter_rows(named=True):
		print(r['junction'] + '\t' + '\t'.join([str(r[c]) for c in cols]))


# For quick-loading into other applications e.g. Circos
def cache_objects():
	d = test_data / 'junction_stats_objects'
	if d.exists():
		print('Object cache already exists')
		return
	
	print('Caching objects')
	d.mkdir()
	pathstats.write_ndjson(d / 'pathstats.jsonl')
	joint_paths.write_ndjson(d / 'joint_paths.jsonl')
	joints_pos.write_ndjson(d / 'joints_pos.jsonl')
	accessory_paths.write_ndjson(d / 'accessory_paths.jsonl')

cache_objects()