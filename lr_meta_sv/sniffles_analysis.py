import polars as pl
from ..util import io, vcfpl
from ..figures import plot_utils as plt

################################################################
# Collect the set of samples & output the calls to sniffles.sh
################################################################

cmd = 'sbatch sniffles_per_patient.sh {} {}.{} ../00.Data/{}_modified.fa reads/{}.fastq.gz'
cmd_self = pl.format(cmd, 'patient', 'sample', 'sample', 'sample', 'barcode')
cmd_ref = pl.format(cmd, 'patient', 'sample', 'ref', 'ref', 'barcode')


samples = (
	pl
	# .read_csv('../data_processing/data/lr_meta_samples.csv') # LR meta samples
	.read_csv('/labs/asbhatt/abehr/s/efm/data_processing/data/lr_meta_samples.csv')
	.filter(pl.col('genome').is_not_null()) # has Efm chromosome
	.select('patient', 'sample', 'barcode')
)
patients = (
	samples
	.group_by('patient').agg('sample') # Group samples per patient
	.filter(pl.col('sample').list.len() > 1) # Remove patients which only have one sample
	.with_columns(ref=pl.col('sample').list.sort().list.first()) # Take the first sample as the "basis"/reference
)

samples = samples.join(patients.drop('sample'), on='patient').sort(by=['patient', 'sample'])

cmds = pl.concat([
	samples.with_columns(cmd=cmd_self)['cmd'],
	samples.with_columns(cmd=cmd_ref)['cmd'],
])

# -> Run Sniffles
io.write_bash_script_cmds(cmds.unique(), 'cmds.sh')

###################################
# Define & import Sniffles output
###################################

with open('../08.FaeciumChromosomes/efm_chromosomes.txt') as f:
	efm_chromosomes = [l.strip() for l in f]

for version in ('self', 'patient'):
	snf = []
	for r in samples.iter_rows(named=True):
		patient = r['patient']
		sample = r['sample']
		ref = r['sample'] if version == 'self' else r['ref']

		for run in ('normal', 'mosaic'):
			name = f"{sample}.{ref}.sniffles{'.mosaic' if run == 'mosaic' else ''}.vcf"
			lf = pl.scan_csv(name, separator='\t', comment_prefix='##', infer_schema=False)
			snf.append(lf.with_columns(patient=pl.lit(patient), sample=pl.lit(sample)))
		
	# Filter to SVs against the Efm chromosomes
	snf = pl.concat(snf).filter(pl.col('#CHROM').is_in(efm_chromosomes)).collect()
	
	# Validate Sniffles format & parse additional columns
	vcfpl.validate_format(snf)
	snf = snf.with_columns(vcfpl.parse_format)
	snf = snf.with_columns(vcfpl.parse_info)

	# Drop duplicate SVs
	snf = snf.unique(subset=['sample', 'POS', 'svtype', 'svlen'])
	assert snf.n_unique(subset=['sample', 'ID']) == len(snf)

	snf = (
		snf
		.with_columns(
			svlen=pl # By default, Sniffles deletions have negative length; reverse this.
				.when(svtype='DEL')
				.then(-1*pl.col('svlen'))
				.otherwise('svlen')
		)
		# Ignore SV <100 bp. This keeps parity with PanGraph, and ignores BND type
		.filter(pl.col('svlen') >= 100)
	)

	snf.write_csv(f'sniffles_concat.{version}.tsv', separator='\t')


###############################
# Process "self" Sniffles
# (within-sample heterogeneity)
# and make a simple plot of this (Fig. 5B)
###############################

snf = pl.read_csv('data/sniffles_concat.self.tsv', separator='\t')

snf['svtype'].value_counts(sort=True)
'''──────┬───────┐
│ svtype ┆ count │
│ ---    ┆ ---   │
│ str    ┆ u32   │
╞════════╪═══════╡
│ INS    ┆ 30    │
│ DUP    ┆ 30    │
│ DEL    ┆ 12    │
│ INV    ┆ 4     │
└────────┴─────'''

import plotly.express as px
pdf = (
	snf
	.group_by('sample').len()
	.join(samples.select('sample'), how='right', on='sample')
	.fill_null(0)
)
sample_order = pdf['sample'].sort().to_list()
fig = px.bar(
	pdf, x='sample', y='len',
	category_orders=dict(sample=sample_order),
	template='plotly_white'
)
plt.common_figure_style(fig, publish=True)
fig.write_image('figs/self.svg')


###############################
# Extract IS insertions over time from 'patient' SVs
###############################

snf = (
	pl
	.read_csv('sniffles_concat.patient.tsv', separator='\t')
	.filter(
		# Only consider insertions/deletions, and look in the size range of 1-2 ISEs
		pl.col('svlen').is_between(750, 3000),
		pl.col('svtype').is_in(['INS', 'DEL'])
	)
	.with_columns(seqid=pl.format('{}_{}', 'sample', 'ID'))
)

with open('patient_indels.fa', 'w') as f:
	for sv in snf.iter_rows(named=True):
		seqid = sv['seqid']
		seq = sv['ALT'] if sv['svtype'] == 'INS' else sv['REF']
		f.write(f'>{seqid}\n{seq}\n')

# Now can run ISEScan/Bakta on the output and process it

