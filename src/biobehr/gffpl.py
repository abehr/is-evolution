import polars as pl
from pathlib import Path
from pyfaidx import Fasta


gff_cols = dict(
	seqid = pl.Utf8,
	source = pl.Utf8,
	type = pl.Utf8,
	start = pl.Int64,
	end = pl.Int64,
	score = pl.Float64,
	strand = pl.Utf8, # +, -, ? (relevant but undetermined), . (irrelevant)
	phase = pl.Int64, # '.' if N/A (for non-CDS features)
	attributes = pl.Utf8
)
col_names_with_sample = ['sample'] + list(gff_cols.keys())

def validate_gff_file(f):
	# Convert to Path if not already
	if isinstance(f, str): f = Path(f)

	# See if an HTS-converted version already exists.
	# TODO: note that my GFFPL parser isn't as strict as HTSeq.GFF_Reader, but still,
	# it can't handle the ##FASTA section that Bakta produces. 
	hts_version = f.parent / f'{f.stem}_htseq_compatible{f.suffix}'
	if not f.stem.endswith('_htseq_compatible') and hts_version.is_file():
		f = hts_version
	
	assert f.is_file(), f'GFF file {f} does not exist!'

	return f

class GFFPL:
	def __init__(self, gff_file_or_files_or_records, samples=None):
		if isinstance(gff_file_or_files_or_records, list):
			self.multi = True
			if samples is not None:
				# List of gff file paths, and list of corresponding sample names (must be the same length)
				gff_files = gff_file_or_files_or_records
				self.lf = pl.concat([GFFPL.parse(f, s) for f,s in zip(gff_files, samples)], how='vertical')
			elif isinstance(gff_file_or_files_or_records[0], dict):
				# List of dicts, where each item contains a "gff" (Path or string path to gff file) and a "sample" (sample name).
				# e.g. if you have tabular format data with gff paths & sample names (other fields are ok and will be ignored).
				records = gff_file_or_files_or_records
				self.lf = pl.concat([GFFPL.parse(r['gff'], r['sample']) for r in records], how='vertical')
			else:
				# Assume it's just a list of gff files, and we will name sample based on filename
				records = [Path(r) for r in gff_file_or_files_or_records]
				self.lf = pl.concat([GFFPL.parse(r, r.stem) for r in records])
		elif isinstance(gff_file_or_files_or_records, pl.DataFrame):
			# DataFrame, with required columns "gff" (string path to gff file) and "sample" (sample name
			assert samples is None
			self.multi = True
			data = gff_file_or_files_or_records.select(['gff', 'sample'])
			self.lf = pl.concat([GFFPL.parse(x[0], x[1]) for x in data.iter_rows()])
		elif isinstance(gff_file_or_files_or_records, dict):
			# Single GFF: dict with "gff" (Path or string path to gff file) and "sample" (sample name).
			assert samples is None, 'Sample should be dict attr if dict record provided'
			gff_file = gff_file_or_files_or_records['gff']
			sample = gff_file_or_files_or_records['sample']
			self.lf = GFFPL.parse(gff_file, sample)
		else:
			# Single GFF: Path or string path to gff file and (optional) sample name.
			assert samples is None or isinstance(samples, str), 'Cannot provide sample list for one gff file'
			gff_file = gff_file_or_files_or_records
			self.lf = GFFPL.parse(gff_file, samples)

		self.lf = self.lf.with_columns(
			ID=pl.col('attributes').str.extract(r'ID=([^;]+)', 1),
			product=pl.col('attributes').str.extract(r'product=([^;]+)', 1)
		)
	
	# Note that these filters are permanent, so returning self is kinda misleading.
	# OTOH it allows you to just build your gffpl nicely. 
	# g = GFFPL(file).contig('contig_3').valid_cds().lf

	def contig(self, seqid):
		self.lf = self.lf.filter(pl.col('seqid') == seqid)
		return self
	
	def valid_cds(self):
		self.lf = self.lf.filter(
			(pl.col('type') == 'CDS') & 
			(pl.col('strand').is_in(['+', '-'])) & 
			pl.col('product').is_not_null()
		)
			
		return self
	
	def cds_rna(self):
		self.lf = self.lf.filter(pl.col('type').is_in(['CDS', 'tRNA', 'rRNA']))
		return self

	def parse(gff_file, sample=None):
		lf = pl.scan_csv(
			validate_gff_file(gff_file), 
			separator='\t', 
			comment_prefix='#', 
			has_header=False, 
			schema=gff_cols, 
			null_values='.'
		)
		if sample is not None:
			lf = lf.with_columns(sample=pl.lit(sample)) # add Sample col
			lf = lf.select(col_names_with_sample) # reorder cols
			# Note we actually may not want to reorder cols so that the initial set of columns is always the same, regardless of whether Sample is included or not. 
			# otherwise, sample should be absolutely required. 
		return lf


# case-insensitive pl.Expr
contains = lambda col_name, substr: pl.col(col_name).str.contains(f'(?i){substr}')
contains_any = lambda col_name, substr_list: pl.col(col_name).str.contains_any(substr_list, ascii_case_insensitive=True)

# Inspect row
def look(df, row=0):
	for x in list(df[row]): print(f'{x.name}:\t{x[0]}')


# Takes a pyfaidx.Fasta object and a pl.df gff row
# NOTE: we shouldn't have to check for valid_cds because we can just ensure that the entire df was pre-filtered to only be valid cds.
# feat should be a GFFPL DataFrame row (i.e. a *named* tuple from doing an iter_rows(named=True)).
# Wrapper for below orf function to make calling easier
def get_orf(genome, feat):
	return slice(genome, feat['seqid'], feat['start'], feat['end'], feat['strand'])


def get_flank(genome, feat, length, left=False, right=False):
	assert left or right, 'left/right flank must be chosen'
	s, e, strand, seqid = int(feat['start']), int(feat['end']), feat['strand'], feat['seqid']
	f = lambda st, en, strnd: slice(genome, seqid, st, en, strnd)
	l = f(s-length, s-1, '+') if strand == '+' and left  else f(e+1, e+length, '-') if strand == '-' and left  else None
	r = f(e+1, e+length, '+') if strand == '+' and right else f(s-length, s-1, '-') if strand == '-' and right else None
	return (l, r) if (left and right) else l if left else r

	
def slice(genome, seqid, start, end, strand):
	if not isinstance(genome, Fasta):
		raise TypeError(f'{genome} must be pyfaidx.Fasta')

	nt = genome[seqid][start-1:end] # note start-1 (gff 1-indeded) vs feature.iv.start (0-indexed)
	if strand == '-': nt = nt.reverse.complement
	nt = nt.seq

	if end > len(genome[seqid]): 
		tail = genome[seqid][:end - len(genome[seqid])]
		if strand == '-': tail = tail.reverse.complement
		nt += tail.seq
	
	# Let's decide if we do or don't want to translate...for now we won't. 
	return nt