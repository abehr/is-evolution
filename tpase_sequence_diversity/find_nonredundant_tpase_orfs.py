from pyfaidx import Fasta
import polars as pl
from collections import defaultdict
import itertools
from ..util import io
from ..util.gffpl import GFFPL, get_orf, orf
from ..util.translate import translate
from ..util.efmpl import contains, efm_tpase_family_names



data = pl.read_csv('samples.csv', null_values='') # Generated from 01_collect_data.py
# TODO: Let's leave out the lr_meta collection for now and come back to it in a bit.
iso = data.filter(pl.col('contig').is_null())
gff = GFFPL(iso).valid_cds().lf

tpase = gff.filter(contains('product', 'transposase'))\
	.with_columns(fam=efm_tpase_family_names)\
	.filter(pl.col('fam').is_not_null())

# Select only the columns we need and join with the genomic data
tpase = tpase.select(['sample','seqid','start','end','strand','product','fam'])
tpase = tpase.collect() # honestly at this point it's not that huge. 
tpase = tpase.join(iso.select(['sample','genome']), on='sample')



nt_counter = itertools.count()
aa_counter = itertools.count()
all_unique_nt = defaultdict(lambda: next(nt_counter))
all_unique_aa = defaultdict(lambda: next(aa_counter))

results = []
i = 0
n = tpase.select('sample').n_unique()
for genome,group in tpase.group_by('genome'):
	genome = Fasta(genome[0])
	for feat in group.iter_rows(named=True):
		orf = get_orf(genome, feat['seqid'], feat['start'], feat['end'], feat['strand'])
		if len(orf) < 500: continue # ignore anything that's way too short
		nt_seqid = all_unique_nt[orf]
		aa_seqid = all_unique_aa[translate(orf)]
		feat['nt_seqid'] = nt_seqid
		feat['aa_seqid'] = aa_seqid
		results.append(feat)
	
	i += 1
	io.print_progress(i, n, cadence=1)

df = pl.DataFrame(results)

# Write information to files
df.write_csv('tpase_per_sample.csv')
with open('all_unique_tpase_nt.fna', 'w') as f:
	for seq,seqid in all_unique_nt.items(): f.write(f'>{seqid}\n{seq}\n')
with open('all_unique_tpase_aa.faa', 'w') as f:
	for seq,seqid in all_unique_aa.items(): f.write(f'>{seqid}\n{seq}\n')