#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// include {genomad} from '../nf_eskape_lr/main.nf'
include {run_isescan} from '../nf_eskape_lr/main.nf'

def help() {
	log.info """
	Example:
	nf run \$s/efm/isl3_expansion_timeline \
		-profile scg -c \$nf/main.config -c \$nf/notify/config \
		-params-file example.yml [--test] [-resume]
	"""
}

if (params.help) {
	help()
	exit 0
}


/*  
Preprocess short read (Illumina) data using less intermediate space
by running multiple things and cleaning up as you go.
Note: I actually seem to find that prefetch & fasterq-dump are FASTER when you do not try to parallelize them.
Or maybe it's the memory allocation amount. Either way, while trimming might be slowed down by this,
it seems overall faster to do it with a single core, plus it's not very slow overall at scale. 
*/
process prefetch {
	tag "${srr}_prefetch"
	cpus 1
	memory '6 GB'
	time '1h'
	errorStrategy 'ignore'

	input:
	val srr

	output:
	tuple val(srr), path("sra/${srr}/${srr}.sra")

	"""
	echo "aab-preproc-prefetch"
	run_prefetch() {
		${params.sra_tools_bin}/prefetch $srr -O sra
		echo "Original prefetch output status: \$?"
		if [ -f sra/${srr}/${srr}.sra ]; then
			echo "Output file created successfully"
			exit 0
		else
			echo "Output file missing"
			exit 1
		fi
	}
	run_prefetch
	"""
}

process preprocess_illumina_paired {
	tag "${srr}_preprocess"
	cpus 1
	memory '6 GB'
	time '1h'
	errorStrategy 'ignore'
	conda '/labs/asbhatt/abehr/envs/dedup_trim'

	input:
	tuple val(srr), path(sra)

	output:
	tuple val(srr), path("${srr}_dedup_R1_val_1.fq.gz"), path("${srr}_dedup_R2_val_2.fq.gz")
	
	"""
	echo "aab-preproc-fasterq-dump"
	${params.sra_tools_bin}/fasterq-dump --split-files $sra
	
	echo "aab-preproc-dedup"
	hts_SuperDeduper -1 ${srr}_1.fastq -2 ${srr}_2.fastq -f ${srr}_dedup
	rm ${srr}_1.fastq ${srr}_2.fastq

	echo "aab-preproc-trim"
	trim_galore --fastqc --paired ${srr}_dedup_R1.fastq.gz ${srr}_dedup_R2.fastq.gz
	rm ${srr}_dedup_R1.fastq.gz ${srr}_dedup_R2.fastq.gz

	echo "aab-preproc-done"
	"""
}

process spades {
	tag "${srr}_spades"
	cpus 8
	memory '32 GB'
	time '1h'
	errorStrategy 'ignore'
	// need publishDir and conda
	publishDir "$params.outdir/${srr}", mode: 'link', saveAs: { it.split('/')[-1] }

	input:
	tuple val(srr), path(r1), path(r2)

	output:
	tuple val(srr), path('spades/contigs.fasta'), emit: assembly
	path "spades/*.log"

	"""
	${params.spades_bin}/spades.py --isolate -t $task.cpus -m 100 \
		-1 ${r1} -2 ${r2} -o spades
	"""
}

process quast {
	tag "${srr}_quast"
	cpus 1
	memory '2 GB'
	errorStrategy 'ignore'
	conda '/labs/asbhatt/abehr/envs/quast'
	publishDir "$params.outdir/${srr}", mode: 'link', saveAs: { it.split('/')[-1] }

	input:
	tuple val(srr), path(assembly)

	output:
	path 'quast_results/latest/report.tsv'

	// could add add'l like --gene-finding --conserved-genes-finding --rna-finding
	"""
	quast.py $assembly -r $params.ref_genome -g $params.ref_annotation --fast
	"""
}

workflow {
	ch_sra_runs = Channel
		.fromPath(params.sra_run_accessions, checkIfExists:true)
		.splitText()
		.map { line -> line.trim() } // trim leading/trailing whitespace
		.filter { it }  // Exclude any empty or whitespace-only lines

	if (params.test) {
		ch_sra_runs = ch_sra_runs.take(2)
	}

	sra = prefetch(ch_sra_runs)
	fastq = preprocess_illumina_paired(sra)
	assemblies = spades(fastq).assembly
	quast(assemblies)
	run_isescan(assemblies)
	// genomad(assemblies)
}