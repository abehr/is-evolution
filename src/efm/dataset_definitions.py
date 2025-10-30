import polars as pl

NA = pl.lit('undefined')

# collections (mirrored by locate)
NCBI = 'ncbi' # all NCBI collections
NCBI_NANOPORE = 'ncbi_nanopore'
NCBI_COMPLETE = 'ncbi_complete' #_genus

SHC = 'shc'
SHC_ISO = 'shc_nanopore' #_set
SHC_ENTE = 'shc_atypical_enterococcus'

ENTE_REF = 'ente_ref'

LR_META = 'lr_meta'

# External-facing collection names; some logic is different
C_NCBI = 'NCBI'
C_NCBI_COMPLETE = 'NCBI-complete'
C_NCBI_LONGREAD = 'NCBI-longread'
C_SHC = 'SHC'
C_SHC_ISO = 'SHC-iso' # incl. atypical ente
C_SHC_META = 'SHC-meta'
C_HC_ENT = 'HC-Ent'

# Filter to HC-Ent logic
HC_ENT = (
	pl.col('collection').is_in([C_NCBI_COMPLETE, C_NCBI_LONGREAD, C_SHC_ISO]),
	pl.col('genus') == 'Enterococcus'
)

# Convert "legacy" dataset names to the simplified external-facing names which make more intuitive sense
# (e.g. for supplementary data tables)
coll = pl.col('collection')
convert_collection_name = (
	pl
	.when(coll.str.starts_with(NCBI_COMPLETE)).then(pl.lit(C_NCBI_COMPLETE))
	.when(coll.str.starts_with(NCBI_NANOPORE)).then(pl.lit(C_NCBI_LONGREAD))
	.when(coll.str.starts_with(SHC_ISO)).then(pl.lit(C_SHC_ISO))
	.when(coll == SHC_ENTE).then(pl.lit(C_SHC_ISO)) # note ERV 165
	.otherwise(NA)
)

survey_genera = [
	'Acinetobacter',
	'Bordetella',
	'Enterobacter',
	'Enterococcus',
	'Escherichia',
	'Klebsiella',
	'Pseudomonas', 
	'Salmonella',
	'Shigella',
	'Staphylococcus',
	'Streptococcus'
]

# NOTE: Removed E. faecalis from the ESKAPEE list during revisions.

eskape_order = [
	'Enterobacter spp.',
	'S. aureus',
	'K. pneumoniae',
	'A. baumannii',
	'P. aeruginosa',
	'E. coli',
	'E. faecium',
	# 'E. faecalis',
]

eskape_order_df = pl.DataFrame(eskape_order, schema=['eskape']).with_row_index(name='eskape_idx')

eskape_taxa = (
	pl
	.when(genus='Enterococcus', species='faecium').then(pl.lit('E. faecium'))
	.when(genus='Enterococcus', species='faecalis').then(pl.lit('E. faecalis'))
	.when(genus='Staphylococcus', species='aureus').then(pl.lit('S. aureus'))
	.when(genus='Klebsiella', species='pneumoniae').then(pl.lit('K. pneumoniae'))
	.when(genus='Acinetobacter', species='baumannii').then(pl.lit('A. baumannii'))
	.when(genus='Pseudomonas', species='aeruginosa').then(pl.lit('P. aeruginosa'))
	.when(genus='Enterobacter').then(pl.lit('Enterobacter spp.'))
	.when(genus='Escherichia', species='coli').then(pl.lit('E. coli'))
	.otherwise(pl.lit(None))
)

eskape_shorthand = (
	pl
	.when(eskape='E. faecium').then(pl.lit('Fm'))
	.when(eskape='E. faecalis').then(pl.lit('Ef'))
	.when(eskape='S. aureus').then(pl.lit('Sa'))
	.when(eskape='K. pneumoniae').then(pl.lit('Kp'))
	.when(eskape='A. baumannii').then(pl.lit('Ab'))
	.when(eskape='P. aeruginosa').then(pl.lit('Pa'))
	.when(eskape='Enterobacter spp.').then(pl.lit('En'))
	.when(eskape='E. coli').then(pl.lit('Ec'))
	.otherwise(pl.lit(None))
)