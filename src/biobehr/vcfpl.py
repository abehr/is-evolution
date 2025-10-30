import polars as pl
from pathlib import Path


# TODO: implement more strict schema for faster i/o, scanning, etc.

def parse_sniffles_vcf(vcf_file, sample=None, ignore_validation=False):
	if isinstance(vcf_file, str): vcf_file = Path(vcf_file)
	validate_file(vcf_file, ignore_validation)
	
	snf = pl.read_csv(vcf_file, separator='\t', comment_prefix='##')
	if len(snf) == 0:
		print(f'Empty VCF file {vcf_file}')
		return pl.DataFrame()
	
	validate_format(snf)

	# Add sample column, if sample is provided, and parse/expand other columns
	snf = snf.with_columns(parse_format)
	snf = snf.with_columns(parse_info)
	snf = snf.drop('FORMAT', 'SAMPLE')
	if sample: snf = snf.with_columns(sample=pl.lit(sample))

	return snf



vcf_file_format = 'VCFv4.2'
sniffles_version = 'Sniffles2_2.6.0'
def validate_file(vcf_file, ignore_validation=False):
	assert vcf_file.is_file(), 'VCF file does not exist'
	assert vcf_file.suffix == '.vcf', 'Not a VCF file'
	
	with open(vcf_file) as f:
		file_format = f.readline().split('=')[1].strip()
		version = f.readline().split('=')[1].strip()
	
	error = ''
	if file_format != vcf_file_format:
		error += f'Expected {vcf_file_format} (found {file_format})\n'
	if version != sniffles_version:
		error += f'Expected {sniffles_version} (found {version})\n'
	
	if not ignore_validation: assert len(error) == 0, error
	
	# If we ignore validation, we still print the error. 
	if len(error) > 0: print(error)

def validate_format(snf):
	assert snf['FORMAT'].n_unique() == 1, 'Multiple formats?'
	assert snf['FORMAT'].unique().item() == 'GT:GQ:DR:DV'


fmt = pl.col('SAMPLE').str.split(':').list
parse_format = [
	fmt.get(0).alias('GT'), # Genotype
	fmt.get(1).cast(pl.Int64).alias('GQ'), # Genotype quality
	fmt.get(2).cast(pl.Int64).alias('DR'), # Number of reference reads
	fmt.get(3).cast(pl.Int64).alias('DV') # Number of variant reads
]



# Info field description for Sniffles: https://github.com/fritzsedlazeck/Sniffles/wiki/Output#info-field-description
# NOTE: not accounting for e.g. mosaic and other types of patterns, not tested on multi results
extract_kv = lambda key, dtype: pl.col('INFO').str.extract(rf'{key}=([^;]+)', 1).cast(dtype).alias(key.lower())
desired_info_keys = [('SVLEN', pl.Int64), ('SVTYPE', pl.String), ('SUPPORT', pl.Int64), ('VAF', pl.Float64)]
parse_info = [extract_kv(key, dtype) for key, dtype in desired_info_keys]