#! /usr/bin/env nextflow

process isescan {
	tag "${id}_isescan"
	cpus 4
	memory '16 GB'
	publishDir "$params.outdir/isescan", mode: 'link', saveAs: { fname ->
		fname.endsWith('.csv') ? "${id}.csv" :
		fname.endsWith('.sum') ? "${id}.sum" :
		fname
	}

	input:
	tuple val(id), path(assembly, name: 'contigs.fasta')

	output:
	path 'isescan/contigs.fasta.csv'
	path 'isescan/contigs.fasta.sum'

	"""
	isescan.py --output isescan --seqfile $assembly --nthread $task.cpus
	"""
}

// Compare all data

process mash_dist_reads {
	tag "${id}_mash_dist"
	memory '2 GB'
	publishDir "$params.outdir/mash_dist", mode: 'link', saveAs: { "${id}.txt" }

	input:
	tuple val(id), path(reads_sketch)
	path ref_sketch

	output:
	path 'mash_dist.txt'

	// Note the -m 2 flag. Wouldn't want to do this for contig/assemblies. 
	"""
	mash dist -m 2 $ref_sketch $reads_sketch > mash_dist.txt
	"""
}


process fastANI {
	tag "${id}_fastANI"
	memory '2 GB'
	publishDir "$params.outdir/fastANI", mode: 'link', saveAs: { "${id}.txt" }

	input:
	tuple val(id), path(assembly)
	path ref_genome

	output:
	path 'fastani.out'

	"""
	fastANI -q $assembly -r $ref_genome -o fastani.out
	"""
}

process checkM {
	tag "${id}_checkM2"
	cpus 4
	memory '45 GB' // this might take some trial & error, can also do --lowmem
	publishDir "$params.outdir/checkM2", mode: 'link', saveAs: { "${id}.txt" }

	input:
	tuple val(id), path(assembly)

	output:
	path 'checkm2/quality_report.tsv'

	"""
	checkm2 predict --threads 4 -i $assembly -o checkm2 --remove_intermediates
	"""
}


workflow {
	ch_assembly = Channel
		.fromPath(params.assembly)
		.map { f -> [f.baseName, f] }
	
	ch_sketch = Channel
		.fromPath(params.mash_sketch)
		.map { f -> [f.baseName, f] }

	// Remember that these need to be "value" channels (not fromPath) so they can be consumed multiple times.
	ch_ref_genome = Channel.value(file(params.ref_genome))
	ch_ref_sketch = Channel.value(file(params.ref_sketch))

	isescan(ch_assembly)
	fastANI(ch_assembly, ch_ref_genome)
	checkM(ch_assembly)
	mash_dist(ch_sketch, ch_ref_sketch)
}
