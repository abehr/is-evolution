import polars as pl
import json
from pyfaidx import Fasta
from tqdm import tqdm
from pathlib import Path
from collections import OrderedDict
import plotly.graph_objects as go

from biobehr.gff_utils import add_ise_family_column


def parse_joints(joints_pos_json):
	joint_data_fields = ['jcb','jab','jae','jce','junction_strand']
	
	return (
		pl
		.read_json(joints_pos_json)
		.unpivot(variable_name='junction')
		.unnest('value')
		.unpivot(index='junction', variable_name='sample')
		.with_columns(pl.col('value').list.to_struct(fields=joint_data_fields))
		.unnest('value')
	)

# gene_overlaps_junction = ((pl.col("end") >= pl.col("jab")) & (pl.col("start") <= pl.col("jae")))
gene_within_junction = (
	(pl.col('start') >= pl.col('jab')) & 
	(pl.col('end') >= pl.col('jae'))
) 
gene_within_junction_spanning_ori = (
	(pl.col('jab') > pl.col('jae')) &
	(
		(pl.col('start') >= pl.col('jab')) | 
		(pl.col('dnd') <= pl.col('jae'))
	)
)


def local_pan_genomes(gff_lf, joints):
	gff_record_types = ['CDS', 'tRNA', 'rRNA']
	return (
		joints#.lazy()
		.join(gff_lf.filter(pl.col('type').is_in(gff_record_types)), on='sample', how='inner')
		.filter(gene_within_junction | gene_within_junction_spanning_ori)
		.group_by(['junction', 'sample'])
		.agg(
			genes=pl.struct(['start', 'end', 'strand', 'product']),
			products=pl.col('product'),
			n_genes=pl.len(),
			junction_strand=pl.col('junction_strand').first() # doesn't matter
		)
		.with_columns(
			pl
			.when(junction_strand=0)
			.then(pl.col('products').list.reverse())
			.otherwise('products')
			.alias('products')
		)
	)#.collect()


def parse_subgraph_paths(subgraphs_dir):
	joint_paths = []
	for subgraph_json in Path(subgraphs_dir).glob('*.json'):
		with open(subgraph_json) as f:
			subgraph_paths = json.loads(f.readline().strip())['paths']

		for path in subgraph_paths:
			joint_paths.append(dict(
				junction = subgraph_json.stem,
				sample = path['name'],
				path = [x['id'] for x in path['blocks']],
				position = path['position']
			))
	
	return pl.DataFrame(joint_paths).with_columns(
		path_str=pl.col('path').list.join('_'),
		path_length=pl.col('path').list.len(),
		non_empty=(pl.col('path').list.len() > 2) # note that path length 1 is ambiguously empty bc it's not a true junction, but those get filtered out anyway. 
	)


def joints_pos_with_spacers(joints_pos, genomes_dir):
	spacer = 'N'*2000
	spacers = []
	for sample,group in joints_pos.group_by('sample'):
		genome = Fasta(genomes_dir / f'{sample[0]}.fa')
		assert len(genome.records) == 1
		record = genome[list(genome.keys())[0]]
		for r in group.iter_rows(named=True):
			if r['jcb'] != None:
				junction_seq = subseq(record, r['jcb'], r['jce'], coords=True)
				accessory_seq = subseq(record, r['jab'], r['jae'], coords=True)
				r['jct_contains_N'] = 'N' in junction_seq
				r['n_spacers_jct'] = junction_seq.count(spacer)
				r['n_spacers_acc'] = accessory_seq.count(spacer)
			spacers.append(r)
	
	return pl.DataFrame(spacers)

# gff is a gffpl df; strand is bool.
def transform_gff(gff, region_start, region_end, junction_strand):
	"""
	Filter and transform the GFF records so that:
	1. Only records fully within [region_start..region_end] (with possible wrap-around)
	   are retained (no partial overlaps).
	2. Coordinates are shifted so that region_start or region_end maps to 1, depending on strand.
	3. The `strand` column is flipped if junction_strand == False.
	"""

	# Convenience expressions
	gene_crosses_origin = (pl.col("start") > pl.col("end"))
	gene_overlaps_region = (
		(~gene_crosses_origin) &
		(pl.col("start") >= region_start) &
		(pl.col("end") <= region_end)
	)
	gene_overlaps_region_crossing_origin = (
		(
			# for genes that don't cross the origin: start >= region_start OR end <= region_end
			(~gene_crosses_origin) &
			((pl.col("start") >= region_start) | (pl.col("end") <= region_end))
		)
		|
		(
			# For genes that do cross the origin: start >= region_start AND end <= region_end
			(gene_crosses_origin) &
			(pl.col("start") >= region_start) & (pl.col("end") <= region_end)
		)
	)

	gff_filtered = gff.filter(gene_overlaps_region) if region_end > region_start else gff.filter(gene_overlaps_region_crossing_origin)

	# ========== SHIFT COORDS ========== 
	if junction_strand:
		# Positive-strand slice
		# new_start = start - region_start + 1
		# new_end   = end   - region_start + 1
		gff_shifted = gff_filtered.with_columns([
			(pl.col("start") - region_start + 1).alias("new_start"),
			(pl.col("end")   - region_start + 1).alias("new_end"),
			# strand remains as is on the positive slice
			pl.col("strand").alias("new_strand")
		])
	else:
		# Negative-strand slice
		# region_end corresponds to new coordinate 1
		# new_start = (region_end - end)   + 1
		# new_end   = (region_end - start) + 1
		# Also flip the strand
		gff_shifted = gff_filtered.with_columns([
			(region_end - pl.col("end")   + 1).alias("new_start"),
			(region_end - pl.col("start") + 1).alias("new_end"),
			pl.when(strand='+').then(pl.lit('-')).otherwise(pl.lit('+')).alias('new_strand')
		])

	# ========== Reorder start/end if needed ========== 
	# Because of the flipping logic, it's possible that new_start might come out
	# larger than new_end. Usually in GFF, start <= end, so we can enforce that:
	gff_final = gff_shifted.with_columns([
		pl.when(pl.col("new_start") <= pl.col("new_end"))
		.then(pl.col("new_start"))
		.otherwise(pl.col("new_end"))
		.alias("start"),

		pl.when(pl.col("new_start") <= pl.col("new_end"))
		.then(pl.col("new_end"))
		.otherwise(pl.col("new_start"))
		.alias("end"),

		pl.col("new_strand").alias("strand"),
		# Keep other columns (like product) if you want
		pl.col("product")
	]).drop(["new_start", "new_end", "new_strand"])

	return gff_final.sort(by='start')


def extract_junction_gff(gff_dict, joints_pos):
	print('Extract GFF annotations for every junction across all genomes')
	junction_gff = []
	for j in tqdm(joints_pos.iter_rows(named=True), total=len(joints_pos)):
		if j['jcb'] == None: continue # ignore for null junction
		g = transform_gff(gff_dict[j['sample']], j['jcb'], j['jce'], j['junction_strand'])
		# Let seqid correspond to junction ID (bc we have sliced out each joint into its own seq)
		g = g.with_columns(seqid=pl.lit(j['junction']), sample=pl.lit(j['sample']))
		junction_gff.append(g)

	# Concatenate, and select only the columns we're interested in for now
	cols = ['seqid','sample','source','type','start','end','strand','product']
	return pl.concat(junction_gff).select(cols)


def junction_accessory_gff(joint_paths, junction_gff):
	accessory_gff = (
		joint_paths
		.filter(pl.col('path_length') > 2) # no need to look at empty paths (by definition)
		.join(
			junction_gff.rename({'seqid':'junction'}),
			on = ['junction', 'sample'],
			how = 'left',
			maintain_order = 'left_right' # make sure that the junction_gff gene order is preserved per-sample
		)
		.filter(
			pl.col('start') >= pl.col('position').list.get(1), # gene starts after the end of the left-flanking core block
			pl.col('end') <= pl.col('position').list.get(-2) # gene ends before the start of the right-flanking core block
		)
	)
	accessory_gff = add_ise_family_column(accessory_gff, col_name='isfam', nofam_category_name='none')
	return accessory_gff



# Get sequence range as a string, even if wraps around the origin
# Note that start and end are zero-indexed by how pyfaidx records are indexed.
# By default we'll assume that your start/end are zero-indeed, i.e. pythonic,
# which is how pyfaidx record indexing works.
# If your start/end are "fasta coordinates", they are 1-indexed. use coords=True.
# Copied over from currently-unused seq util script
def subseq(record, start, end, coords=False):
	if coords: start -= 1
	if start < end:
		return record[start:end].seq
	elif start > end:
		part1 = record[start:len(record)].seq # assume circular sequence
		part2 = record[0:end].seq # from ori to end
		return part1 + part2
	else:
		return '' # no sequence if start==end


# Make the correct Sankey objects to be ingested into go.Sankey
# from an easier-to-grok 2-layered dict. 
def sankey(sankey_dict):
	labels_ordered = OrderedDict()
	for src, targets in sankey_dict.items():
		labels_ordered[src] = None
		for tgt in targets.keys():
			labels_ordered[tgt] = None

	labels = list(labels_ordered.keys())
	label2idx = {label: i for i, label in enumerate(labels)}

	# 2) Build source, target, and value lists.
	sources = []
	targets = []
	values = []
	for src, subdict in sankey_dict.items():
		for tgt, val in subdict.items():
			sources.append(label2idx[src])
			targets.append(label2idx[tgt])
			values.append(val)

	fig = go.Figure(data=[go.Sankey(
		arrangement='freeform',
		node=dict(label=labels),
		link=dict(source=sources, target=targets, value=values)
	)])

	return fig
