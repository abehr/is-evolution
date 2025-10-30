import sys
import subprocess
from pathlib import Path

def preprocess_gff_for_htseq_and_gffpl(input_gff, output_gff=None, force=True):
	"""
	HTSeq has a much more performant GFF parser (vs. BCBio, gffutils, gffpandas, etc.)
	However, it can't handle when strand is '?', even though this is valid GFF3. 
	We also trim off the fasta part of the GFF, if it exists, because it's redundant to save it again.
	"""
	i = Path(input_gff)
	assert i.is_file(), f'Input GFF {i} does not exist!'
	assert i.suffix in ('.gff', '.gff3'), f'Input file {i} does not look like GFF!'

	if output_gff is None:
		output_gff = i.parent / f'{i.stem}_htseq_compatible{i.suffix}'
	else:
		output_gff = Path(output_gff)
	
	if not force:
		assert not output_gff.exists(), f'Output GFF {output_gff} already exists!'

	awk_script = '''
	BEGIN { FS="\t"; OFS="\t"; }
	/^##FASTA/ { exit }
	/^#/ { print; next }
	NF != 9 { print; next }
	$7 == "?" { $7 = "." }
	{ print }
	'''

	try:
		cmd = ['awk', awk_script, str(i)]
		with open(output_gff, 'w') as outfile:
			subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
		return output_gff
	except subprocess.CalledProcessError as e:
		print(f"Error in processing GFF file: {e.stderr}")
		raise e
	
if __name__ == '__main__':
	input_gff = sys.argv[1]
	output_gff = sys.argv[2] if len(sys.argv) > 2 else None
	preprocess_gff_for_htseq_and_gffpl(input_gff, output_gff)