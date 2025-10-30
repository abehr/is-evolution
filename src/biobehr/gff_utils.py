import polars as pl
import re

tpase_annotation = pl.col('product').str.contains('(?i)transposase')

# This is a heuristic, has been tested mostly on Bakta-style annotations,
# won't work perfectly. Should always manually inspect outputs. Almost always,
# Bakta annotation prints the IS family name first, then a more specific name.
# PGAP may need additional attention to properly escape terms.
extract_tpase_fam = (
	pl.col('product') 
		.str.replace_all(r"(-|%2C|[Tt]ransposase)", ' ')
		.str.strip_chars()
		.str.split(' ')
		.list.first()
)

# Can work on a DataFrame or LazyFrame GFF
def filter_IS_tpase(gff):
	return (
		gff
		.filter(tpase_annotation) # tpase annotations only
		.with_columns(fam=extract_tpase_fam) # find IS family name, if possible
		.filter(pl.col('fam').str.starts_with('IS')) # filter to likely known IS annotations
	)

# Same as above, but just add the column and don't filter, and a couple other params that you can tweak
# (primarily for compatibility with the Pangraph utils).
def add_ise_family_column(gff, col_name='fam', nofam_category_name=''):
	return (
		gff
		# Assign tpase annotations to their IS family
		.with_columns(pl.when(tpase_annotation).then(extract_tpase_fam).otherwise(pl.lit(None)).alias(col_name))
		# Assign a category for tpase annotations that are not specifically IS tpases
		.with_columns(
			pl
			.when((pl.col(col_name) == '') | ~pl.col(col_name).str.starts_with('IS'))
			.then(pl.lit(nofam_category_name))
			.otherwise(col_name)
			.alias(col_name)
		)
	)


'''
Extract IS name from the annotation more specifically than the family level, if possible.
For example:
- "IS256 family ISEf1 transposase" -> "ISEf1"
- "IS30 family transposase" -> ""

Input a df, and a list of IS family names that are present in the annotations (to mask).
Finds the first word in the annotation that starts with "IS" that isn't in the set of fams.
'''
def extract_IS_name_from_annotation(family_names):
	# Reverse sort fam list by length (s.t. e.g. "IS30" and "IS3" are handled correctly)
	fams_escaped_sorted = list(map(re.escape, sorted(family_names, key=lambda x: -len(x))))
	fams_pattern = r'\b(?:' + '|'.join(fams_escaped_sorted) + r')\b' # regex word boundaries
	return (
		pl.col('product')
		.str.replace_all(r"(-|%2C|[Tt]ransposase)", ' ') # remove as above
		.str.replace_all(fams_pattern, ' ') # remove the IS family name
		.str.strip_chars().str.replace_all(r"\s+", "_") # replace any-length whitespace with underscore
		.str.split('_') # then split on underscore
		.list.eval(pl.element().filter(pl.element().str.contains('^IS'))) # words that start with "IS"
		.list.first().fill_null('') # take the first one, if it exists.
	)