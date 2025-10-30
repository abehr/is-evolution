import sys
import polars as pl
from pathlib import Path
from collections import defaultdict
import itertools
from pyfaidx import Fasta
import plotly.express as px
import math
from tqdm import tqdm

from efm import dataset_definitions as defs, plot_utils as plt
from efm.config import cfg, data
from biobehr import gff_utils, translate
from biobehr.gffpl import GFFPL, get_orf, get_flank

def main():
	samples, tpase_gff = get_sample_list_and_all_tpase_annotations()
	all_unique_tpase_seqs, stats_per_tpase_variant = extract_all_tpase_sequences_from_annotations(samples, tpase_gff)
	fig = cumulative_distribution_log_plot(stats_per_tpase_variant)
	fig.write_image(data.output / '01_tpase_copies.svg')

	# Find the set of most common tpase seqs (optional: write to fasta & CSV)
	most_common_tpase = get_most_common_tpase_set(all_unique_tpase_seqs, stats_per_tpase_variant)

	# Find flanks and compute boundaries for these (optional: write to CSV)
	extract_tpase_flanking_sequence_context(samples, tpase_gff, most_common_tpase)


def get_sample_list_and_all_tpase_annotations():
	if not cfg.source.genome.is_dir() or not cfg.source.gff.is_dir():
		sys.exit('Missing required config dirs: source.genome & source.gff')

	nf = len([x for x in Path(cfg.source.genome).glob('*.fna')])
	ng = len([x for x in Path(cfg.source.gff).glob('*.gff3')])

	if nf == 0 or ng == 0:
		sys.exit('Downloaded genome fasta & GFFs (NCBI and SHC) are required for this analysis')

	print(f'Get sample list & all transposase annotations (source dirs contain {nf} fasta & {ng} gff files)')
	
	cols = ['sample', 'collection', 'genus', 'species', 'genome_length', 'contigs']
	samples = (
		pl
		# Select isolate genomes
		.concat([
			pl.read_csv(data.ncbi.combined_qc_filtered).select(cols),
			pl.read_csv(data.shc.iso).select(cols)
		])
		# Select E. faecium
		.filter(species='faecium')
		# Only use NCBI 'complete' and SH genomes to avoid potential issues from older nanopore chem,
		# which is ideal for this analysis where exact gene orf nucleotide sequences are compared.
		.filter(
			pl.col('collection').str.starts_with(defs.NCBI_COMPLETE) | 
			pl.col('collection').str.starts_with(defs.SHC_ISO)
		)
		# Add path to genome & bakta annotation files based on location provided in config.
		# This assumes that the files are named by sample.
		.with_columns(
			genome=pl.format('{}/{}.fna', pl.lit(str(cfg.source.genome)), 'sample'),
			gff=pl.format('{}/{}.gff3', pl.lit(str(cfg.source.gff)), 'sample')
		)
	)

	# This consists of 593 samples (493 NCBI complete; 100 SHC).
	# There are only 465 of 493 NCBI complete Efm that are in the "combined final samples", 
	# where ISEScan or geNomad may have failed. But we can include anything here with a genome+gff.
	gff = GFFPL(samples).valid_cds().lf

	tpase_gff = (
		gff_utils.filter_IS_tpase(gff)
		.select('sample','seqid','start','end','strand','product','fam')
		.collect()
	)
	print(f'Found {len(tpase_gff)} IS transposases across {len(samples)} samples')
	return samples, tpase_gff

# Now that we have the list of all IS family tpase annotations, extract them from each genome, and compute their distribution across all the samples
def extract_all_tpase_sequences_from_annotations(samples, tpase, n_most_common=6):
	nt_counter = itertools.count()
	# aa_counter = itertools.count()
	all_unique_nt = defaultdict(lambda: next(nt_counter))
	# all_unique_aa = defaultdict(lambda: next(aa_counter))

	results = []
	# n = tpase.select('sample').n_unique()
	# i = 0
	tpase_per_genome = tpase.join(samples['sample', 'genome'], on='sample', how='left').group_by('genome')
	for genome,group in tqdm(tpase_per_genome, total=tpase.select('sample').n_unique()):
	# for genome,group in tpase.join(samples['sample', 'genome'], on='sample', how='left').group_by('genome'):
		genome = Fasta(genome[0])
		for feat in group.iter_rows(named=True):
			orf = get_orf(genome, feat)
			# Ignore anything that's clearly too short. 600 (~200aa) would already be pretty 
			# safe cutoff, but we'll go down to ~167 to be even more conservative. 
			if len(orf) < 500: continue
			nt_seqid = all_unique_nt[orf]
			# aa_seqid = all_unique_aa[translate(orf)]
			feat['nt_seqid'] = nt_seqid
			# feat['aa_seqid'] = aa_seqid
			results.append(feat)

	df = pl.DataFrame(results).drop('genome')

	nt_seq = pl.DataFrame(list(all_unique_nt.items()), schema=[('seq', str), ('nt_seqid', int)], orient='row')

	# Write information to files (if you want -- not necessary anymore with fast local processing)
	# df.write_csv('efm_tpase_per_sample.csv')
	# with open('efm_all_unique_tpase_nt.fna', 'w') as f:
		# for seq,seqid in all_unique_nt.items(): f.write(f'>{seqid}\n{seq}\n')
	# with open('efm_all_unique_tpase_aa.faa', 'w') as f:
		# for seq,seqid in all_unique_aa.items(): f.write(f'>{seqid}\n{seq}\n')

	# Average tpase count per genome = 124
	# (previously 112, but now we're including more tpase families)
	# df.group_by('sample').len().mean()

	# Subset to the most common IS families to plot (though you could look at more)
	fams = tpase['fam'].value_counts(sort=True)[:n_most_common]['fam'].to_list()

	# Compute the stats per exact IS tpase copy
	stats_per_tpase_variant = (
		df
		.filter(pl.col('fam').is_in(fams))
		.group_by('nt_seqid')
		.agg(
			fam=pl.col('fam').mode().first(),
			product=pl.col('product').mode().first(),
			count=pl.len()
		)
		.with_columns(avg_count=pl.col('count') / df.select('sample').n_unique())
		.with_columns(pdt=gff_utils.extract_IS_name_from_annotation(fams))
		.sort(by='avg_count', descending=True)
	)

	print(f'Found {len(stats_per_tpase_variant)} unique tpase variant across the {n_most_common} most common IS families')
	return nt_seq, stats_per_tpase_variant

# Get the set of common (mean copies per genome >= 1) tpase variants
def get_most_common_tpase_set(all_unique_tpase_seqs, stats_per_tpase_variant, write=True):
	most_common_tpase = (
		stats_per_tpase_variant
		.filter(pl.col('avg_count') >= 1)
		.join(all_unique_tpase_seqs, how='left', on='nt_seqid')
		.with_columns(
			uid=pl
				.when(pl.col('pdt').str.len_chars() > 0)
				.then(pl.format('{}-{}-{}', 'nt_seqid', 'fam', 'pdt'))
				.otherwise(pl.format('{}-{}', 'nt_seqid', 'fam'))
		)
	)

	print(f'Found {len(most_common_tpase)} unique tpase variants with mean count per genome of >= 1')

	if write:
		print('Writing most common tpase seqs to fasta & stats to CSV')
		# Write info out to CSV
		most_common_tpase.drop('seq').write_csv(data.output / 'efm_most_common_tpase_copies.csv')

		# Write sequences out to Fasta
		fasta = (
			most_common_tpase
			.select(pl.format('>{}_n{}\n{}', 'uid', 'count', 'seq'))
			.to_series().str.join('\n').item()
		)
		with open(data.output / 'efm_most_common_tpase_copies.fna', 'w') as f:
			f.write(fasta)
	
	most_common_tpase = most_common_tpase.with_columns(
		aa = pl.col('seq').map_elements(lambda seq: translate.translate(seq), return_dtype=str)
	)
	most_common_tpase.group_by('aa').agg('uid', pl.col('count').sum()).sort(by='count', descending=True)
	# ISEF1 7,36 and IS30 8,57 are duplicated at the amino acid level.
	return most_common_tpase

# Extract the sequence context on either side of every observed tpase copy, for use in determining the flank lengths
def extract_tpase_flanking_sequence_context(samples, tpase_gff, most_common_tpase_variants, write=True):
	print('Extract 500bp of sequence context on either side of every copy of each common tpase')
	mct_nt = {seq:seqid for seq,seqid in most_common_tpase_variants['seq', 'nt_seqid'].iter_rows()}
	nt_counter = itertools.count()
	all_unique_flanks = defaultdict(lambda: next(nt_counter))
	results = []
	tpase_per_genome = tpase_gff.join(samples['sample', 'genome'], on='sample', how='left').group_by('genome')
	for genome,group in tqdm(tpase_per_genome, total=tpase_gff['sample'].n_unique()):
		genome = Fasta(genome[0])
		for feat in group.iter_rows(named=True):
			orf = get_orf(genome, feat)
			if len(orf) < 500: continue # ignore anything that's way too short
			if orf not in mct_nt: continue # ignore anything that's not in the list of 17 most commmon
			
			feat['nt_seqid'] = mct_nt[orf]
			feat['l_flank'], feat['r_flank'] = get_flank(genome, feat, 500, left=True, right=True)
			feat['lfid'], feat['rfid'] = all_unique_flanks[feat['l_flank']], all_unique_flanks[feat['r_flank']]
			results.append(feat)

	df = (
		pl
		.DataFrame(results)
		.drop('genome')
		# Ignore cases where we don't have the entire flank context (there are just a few of those, we lose like 50 out of 36k)
		.filter(pl.col('l_flank').str.len_chars() == 500, pl.col('r_flank').str.len_chars() == 500)
		# Create a unified ID for the combination of left & right flanks, tho actually this doesn't really matter
		# .with_columns(fid=pl.format('{}_{}', 'lfid', 'rfid'))
	)

	# Filter any duplicated flanks; this allows us to do a better motif analysis bc we're not enriching for any more-common positions.
	# This reduces the dataframe size from 36.6k to 5257 (with flank length 100), or 6199 (with flank length 150)
	# and it's 10106 with 200 so it seems fine to do that, we're clearly not getting a huge differnece. 
	flanks = pl.concat([
		df.select('nt_seqid', flank='l_flank', fid='lfid', fpos=pl.lit('L')),
		df.select('nt_seqid', flank='r_flank', fid='rfid', fpos=pl.lit('R'))
	])

	# Filter dupes
	flanks = flanks.group_by(flanks.columns).agg(count=pl.len())
	if write: flanks.write_csv(data.output / 'efm_most_common_tpase_flanks.csv')
	return flanks

# ============ Plot results ====================

# Original plot style: bars, average copy per genome
def original_bar_plot(stats_per_tpase_variant):
	fam_order = stats_per_tpase_variant.group_by('fam').sum().sort(by='count', descending=True)['fam'].to_list()
	fig = px.bar(
		stats_per_tpase_variant, x='avg_count', y='fam', orientation='h',
		category_orders=dict(fam=fam_order),
		hover_data='product',
		template='plotly_white',
		color_discrete_sequence=['#E4C9E1']
	)
	
	fig.update_traces(text=stats_per_tpase_variant['pdt'], textposition='inside')
	fig.update_traces(textfont=dict(size=10, family='Helvetica'))
	fig.update_layout(uniformtext_minsize=10, uniformtext_mode='hide')
	plt.common_figure_style(fig, publish=True)
	fig.update_traces(marker_line_color="black", marker_line_width=0.1)
	
	# Fix the fig dimensions to make sure the labels render how we want, and write out
	fig.update_layout(height=350, width=720)
	return fig

# cumulative distribution (rarefaction curve-ish style) plot
def cumulative_distribution_log_plot(stats_per_tpase_variant):
	pdf = stats_per_tpase_variant.with_columns(
		cumulative=pl.col('count').cum_sum().over('fam'),
		idx=pl.col('count').cum_count().over('fam')
	)

	epsilon = 0.53 # on a log scale, [0.5,1] is the same distance as [1,2] -- slight adjustment bc of marker size & margins etc. 
	zeros = pdf.select('fam').unique().with_columns(cumulative=pl.lit(0), idx=pl.lit(epsilon))
	pdfz = pl.concat([zeros, pdf], how='diagonal_relaxed')
	
	# fig = px.line(pdfz, x='idx', y='cumulative', color='fam', markers=True, color_discrete_map=IS_fam_color_map)
	
	# log version
	fig = px.line(pdfz, x='idx', y='cumulative', color='fam', markers=True, color_discrete_map=plt.IS_fam_color_map)
	fig.update_layout(template='plotly_white', showlegend=True)
	plt.common_figure_style(fig, publish=True)

	fig.update_xaxes(type='log', range=[math.log10(0.5), math.log10(1100)],
		tickvals=[epsilon, 1 , 2 , 5 , 10 , 20 , 50 , 100 , 200 , 500 , 1000], 
		ticktext=[    '0','1','2','5','10','20','50','100','200','500','1000'],
		title_text='No. Unique IS Tpase ORFs' # "sequence copies, nucleotide level"
	)
	fig.update_yaxes(range=[0, 22000], title_text='No. IS Tpase ORFs') # copies, across 600 Efm genomes
	
	fig.update_layout(
		title_text = 'Frequency/distribution of unique IS transposase variants',
		xaxis_title = 'Cumulative no. unique IS transposase variants',
		yaxis_title = 'Cumulative total no. IS transposase copies observed (across 593 E. faecium genomes)'
	)
	return fig



if __name__ == '__main__':
	main()