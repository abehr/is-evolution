#!/usr/bin/env python

import sys
from pathlib import Path
import polars as pl

low_threshold = 20 # Throw out if <20x coverage
high_threshold = 150 # Downsample if > 150x coverage
coverage_target = 110 # Downsample to 100-120x coverage

def main(genome_length_mbp, multiqc_fastqc):
	est_genome_size = float(genome_length_mbp)*1e6
	multiqc_fastqc = Path(multiqc_fastqc)
	assert multiqc_fastqc.is_file()
	assert est_genome_size < 1e7 # Just in case user accidentally provides bp rather than Mbp

	df = (
		pl
		.read_csv(multiqc_fastqc, separator='\t')
		.with_columns(total_length = pl.col('Total Sequences')*pl.col('avg_sequence_length'))
	)

	# Estimate coverage from total (combined R1 & R2) seq length
	est_cov = df['total_length'].sum() / est_genome_size

	# Downsample fraction to get to the target coverage. This number is invalid if coverage is already ok.
	downsample_fraction = coverage_target / est_cov

	# Gate coverage threshold
	cov_status = 'low' if est_cov < low_threshold else 'high' if est_cov > high_threshold else 'ok'

	print(f'{cov_status}\t{downsample_fraction:.2f}')
	return (cov_status, downsample_fraction)










if __name__ == '__main__':
	if len(sys.argv) != 3: sys.exit('Usage: <script> <est_genome_size_mbp> <multiqc_fastqc_txt>')
	main(sys.argv[1], sys.argv[2])
