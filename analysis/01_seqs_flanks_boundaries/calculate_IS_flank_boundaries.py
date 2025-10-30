import sys
import polars as pl
import subprocess
from tqdm import tqdm
import tempfile
from pathlib import Path
import shutil

from efm.config import data
from biobehr import align


'''
This script is run *after* exact_copy_analysis.py (it relies on those outputs).
'''

if shutil.which("mafft") is None:
    sys.exit("Error: 'mafft' not found in PATH. Please install MAFFT or add it to your PATH.")

mct = pl.read_csv(data.output / 'efm_most_common_tpase_copies.csv')
flanks = pl.read_csv(data.output / 'efm_most_common_tpase_flanks.csv')

# Not the prettiest way of handling the temporary files but it's fine. 
tmpdir = Path(tempfile.mkdtemp())

def align_consensus(df, min_length):
	fa_in = tmpdir / 'x.fa' # temp input fasta for aligner
	fa_out = tmpdir / 'a.fa' # temp output aligned fasta

	df_fasta = df.select(pl.format('>{}f_{}n\n{}', 'fid', 'count', 'flank'))
	with open(fa_in, 'w') as f:
		f.write(df_fasta.to_series().str.join('\n').item())

	subprocess.run(
		['mafft', fa_in], 
		stdout=open(fa_out,'w'), stderr=subprocess.PIPE, 
		universal_newlines=True, check=True
	)
	# This kinda depends, would be better to set a higher min region but for is30 it's very small flank. 
	return align.consensus_sequence(fa_out, entropy_threshold=0.5, min_region_size=min_length)


result = []
for tp in tqdm(mct.iter_rows(named=True), total=len(mct)):
	for flank_side in ('L', 'R'):
		c = align_consensus(flanks.filter(nt_seqid=tp['nt_seqid'], fpos=flank_side), min_length=5)
		if c['alignment_length'] == 0:
			print(f"No consensus found for {tp['uid']} {flank_side}")

		c['nt_seqid'] = tp['nt_seqid']
		c['fpos'] = flank_side
		result.append(c)

df = pl.DataFrame(result)
# df.write_csv('consensus.csv')
shutil.rmtree(tmpdir)

# Now, pivot this around a bit so that you can have one row per nt_seqid and join it with the main one. 
consensus = (
	df
	.unpivot(index=['nt_seqid', 'fpos'])
	.with_columns(col_name=pl.format('{}_{}', 'variable', 'fpos')) # give L and R diff column names so it can be pivoted
	.pivot(index='nt_seqid', on='col_name', values='value')
)

consensus = mct.join(df, on='nt_seqid')
consensus.write_csv(data.output / 'consensus_flanks_per_tpase.csv')


