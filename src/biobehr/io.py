import polars as pl
from pathlib import Path
from tqdm import tqdm
import subprocess
import sys
import os
import fnmatch

def validate_file(f, exists=True):
	if isinstance(f, str): f = Path(f)
	if exists: assert f.is_file()
	assert f.suffix in ('.csv', '.tsv', '.txt')
	return f

def validate_outf(f, ext=None):
	if isinstance(f, str): f = Path(f)
	assert not f.exists()
	if ext: assert f.suffix == ext
	return f

# Return the separator to use based on the file extension.
# NOTE: this lambda assumes that f has already been validated by validate_file
sep = lambda f: ',' if f.suffix == '.csv' else '\t'


# Note: for now, only allow whitespace sep
def read_dict(file, header=False):
	file = validate_file(file)
	with open(file) as f:
		if header: f.readline()
		return {x[0]:x[1] for x in (l.strip().split() for l in f)}


'''
Ingest one or multiple input files into a pl.LazyFrame.
Optionally, only ingest a subset of columns, which can also be optionally renamed.
The input file(s) can be csv, tsv, or txt. If txt, it is treated as tsv. 
If multiple input files, they should ofc all have the same column set.

src (input file/files) can be represented as:
- a 1-column df or a series, in which case it is treated as a list of strings corresponding to filepaths
- a list, set, or tuple, in which case it is treated as a list of pathlib.Path objects
- a single pathlib.Path object
- a string, which represents a path to a single file (most simple case)
'''
def pl_read(src, cols=None, dedup=True):
	if isinstance(src, pl.DataFrame):
		# List of paths can also be a pl.DataFrame column or a pl.Series
		assert src.shape[1] == 1, 'Only one column must be provided'
		src = src.to_series()
	if isinstance(src, pl.Series):
		# We will assume, if it's coming in as a pl object, that it is not a pathlib object
		# Note that pl.Series cannot store pathlib.Paths, therefore we can assume it's a series of strings.
		assert src.dtype == pl.String, 'For pl.Series input, expected item dtype to be string'
		src = [Path(x) for x in src]
	if isinstance(src, str):
		# If the input is a string, assume it's a single path to a single file
		src = Path(src)

	if isinstance(src, list) or isinstance(src, tuple) or isinstance(src, set):
		lf = pl.concat([pl.scan_csv(f, separator=sep(f)) for f in src])
	else:
		lf = pl.scan_csv(src, separator=sep(src))
	
	if cols:
		if isinstance(cols, dict):
			lf = lf.select(list(cols.keys())).rename(cols)
		else:
			lf = lf.select(cols)

	return lf.unique() if dedup else lf

def pl_concat(series, fn):
	records = []
	for x in tqdm(series):
		records.append(fn(x))
	
	if isinstance(records[0], dict):
		return pl.DataFrame(records)
	if isinstance(records[0], pl.DataFrame):
		return pl.concat(records)
	
	assert False, f'Cannot concatenate records of type {type(records[0])}'


# stdout, count_lines, and print_progress copied from utils from secreted-peptides repo

def stdout(s):
	sys.stdout.write(s)
	sys.stdout.flush()


def count_lines(fp):
	fp = Path(fp)
	assert fp.is_file(), f'{fp} does not exist.'
	cmd = f'zcat {fp} | wc -l' if fp.suffix == '.gz' else f'wc -l {fp}'
	r = subprocess.run(cmd, capture_output=True, text=True, shell=True)
	assert r.returncode == 0, f'Linecount for {fp} failed ({r.returncode}): {r.stdout} {r.stderr}'
	return int(r.stdout.strip().split()[0])


# Note that os.scandir does not recursively search. So the pattern can only be in the parent dir itself.
def fast_glob(base_dir, pattern, dirs_only=False, files_only=False):
	assert not (dirs_only and files_only), 'Cannot search for only dirs and only files'
	assert '/' not in pattern, 'This method cannot scan subdirectories of the base dir'
	matches = []
	with os.scandir(base_dir) as it:
		for entry in it:
			if fnmatch.fnmatch(entry.name, pattern):
				if dirs_only and not entry.is_dir(): continue
				if files_only and not entry.is_file(): continue
				matches.append(Path(entry.path))
	
	return matches


def write_bash_script_cmds(series, file, overwrite=False):
	if isinstance(file, str): file = Path(file)
	if file.exists() and not overwrite:
		return 'File already exists!'
	
	assert file.suffix == '.sh', 'File is not script!'
	
	cmds = series.str.join('\n').item()
	
	with open(file, 'w') as f:
		f.write('#!/bin/bash\n')
		f.write(cmds + '\n')
	
	os.chmod(file, 0o740) # give me exec permissions


# polars file-column operations (originally in locate module, renamed)

is_file = lambda x: os.path.isfile(x)

'''
Slow (element-wise) filter to check whether files exist at this column.
If a given file doesn't exist, replace with null.
Default column name is "file".
'''
col_exists = lambda col='file': (
	pl
	.when(pl.col(col).map_elements(is_file, return_dtype=pl.Boolean))
	.then(col)
	.otherwise(pl.lit(None))
)

# Check whether files exist at column c, (default column name = "file"),
# *and* filter out the ones that don't exist.
filter_col_files_exist = lambda df, c='file': (
	df
	.with_columns(col_exists(c))
	.filter(pl.col(c).is_not_null())
)

def resolve_files(df, resolve_fp_expr):
	cols = {'sample', 'collection', 'assembly'}
	assert len(cols - set(df.columns)) == 0, 'Input df missing required column'
	# Add 'file' col with path resolved by the input expression;
	# then element-wise check whether each file exists;
	# finally, remove rows in df where the file doesn't exist.
	df = filter_col_files_exist(df.with_columns(file=resolve_fp_expr))

	# Return a series of structs with (sample,file) pairs
	return df.with_columns(data=pl.struct('sample', 'file')).select('data')


'''
Check if column exists. 

e.g., to filter a df down to gff that exist:
df = df.with_columns(gff=locate.gff)
gff_exists = fast_col_exists(df['gff']).filter(exists=True).drop('exists')

df = df.join(gff_exists, how='right', on='gff')

ofc you can do this in other ways. 
'''
from concurrent.futures import ThreadPoolExecutor
def fast_col_exists(df_col: pl.Series, new_col='exists', threads=16):
	paths = df_col.to_list()

	def _exists(p):
		try:
			Path(p).stat()
			return True
		except (FileNotFoundError, OSError):
			return False
	
	with ThreadPoolExecutor(max_workers=threads) as ex:
		exists_flags = list(ex.map(_exists, paths))
	
	return pl.DataFrame([df_col, pl.Series(new_col, exists_flags)])


# Other Polars helper functions

# For column `col`, collapse values not in `common` to "other"
def collapse_other(column, common_values, collapse_to='other'):
	return (
		pl
		.when(pl.col(column).is_in(common_values))
		.then(column)
		.otherwise(pl.lit(collapse_to))
	)