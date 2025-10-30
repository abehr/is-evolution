#! /usr/bin/env nextflow

// nextflow.preview.output = true

/*  
Preprocess short read (Illumina) data using less intermediate space
by running multiple things and cleaning up as you go.
Note: I actually seem to find that prefetch & fasterq-dump are FASTER when you do not try to parallelize them.
Or maybe it's the memory allocation amount. Either way, while trimming might be slowed down by this,
it seems overall faster to do it with a single core, plus it's not very slow overall at scale. 
*/
process prefetch {
	tag "${srr}_prefetch"

	input:
	val srr

	output:
	tuple val(srr), path("sra/${srr}/${srr}.sra")

	"""
	echo "aab-preproc-prefetch"
	run_prefetch() {
		prefetch $srr -O sra
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


/*  Preprocessing
trim_galore: --length default is 20. For paired-end files, both reads of a read-pair
therefore need to be >20. --retain-unpaired allows you to keep single reads that are >20
when their pair is not.
It supports multiple cores, but only for certain parts of the process; most efficient 
across many samples is prob to use single core. 

# Combining unpaired reads:
Even if there are no unpaired, the above will still make an empty file (good for having expected output).
Originally I had zcat *_unpaired_*.fq.gz but this caused a weird loop where the input was re-ingested continually...

# Compress:
pigz multithreaded would be way faster, but this allows us to keep everything in a single process
and remove intermediate outputs while still passing expected final outputs for Nextflow to track. 
*/ 
process preprocess_illumina_paired {
	tag "${srr}_preprocess"
	publishDir "$params.outdir/msh", pattern: "reads.msh", mode: 'link', saveAs: { "${srr}.msh" }

	input:
	tuple val(srr), path(sra)

	output:
		tuple val(srr), path("trim_val_1.fq.gz"), path("trim_val_2.fq.gz"), path("unpaired.fq.gz"), emit: reads
		tuple val(srr), path("trim_val_1_fastqc.zip"), path("trim_val_2_fastqc.zip"), emit: fastqc
		tuple val(srr), path('reads.msh'), emit: msh
		path("stats.log") // dedup stats, if you want them 
	
	"""
	echo "aab-preproc-fasterq-dump"
	fasterq-dump --split-files $sra
	
	echo "aab-preproc-dedup"
	hts_SuperDeduper --uncompressed -1 ${srr}_1.fastq -2 ${srr}_2.fastq -f dedup

	echo "aab-preproc-trim"
	trim_galore --basename trim --paired --fastqc --retain_unpaired dedup_R1.fastq dedup_R2.fastq
	
	echo "aab-preproc-combine-unpaired"
	cat trim_R1_unpaired_1.fq trim_R2_unpaired_2.fq 2>/dev/null > unpaired.fq

	echo "aab-preproc-mash-sketch"
	cat trim_val_1.fq trim_val_2.fq unpaired.fq | mash sketch -m 2 -o reads.msh -

	echo "aab-preproc-compress"
	gzip trim_val_1.fq
	gzip trim_val_2.fq
	gzip unpaired.fq
	
	echo "aab-preproc-cleanup"
	rm *.fastq *.fq

	echo "aab-preproc-done"
	"""
}


// Actually -- I don't think we need the "_" split, bc the result is just filename. 
process multiqc {
	tag "${srr}_multiqc"
	publishDir "$params.outdir/multiQC", mode: 'link', saveAs: { "${srr}_${it.split('/')[-1].split('_')[-1]}" }

	input:
	tuple val(srr), path(z1), path(z2)

	output:
	tuple val(srr), path('multiqc_report.html'), path('multiqc_data/multiqc_data.json'), path('multiqc_data/multiqc_fastqc.txt')

	"""
	multiqc *_fastqc.zip
	"""
}


process eval_coverage {
  tag "${srr}_coverage_eval"
  label 'low_mem'
  publishDir "$params.outdir/coverage_filt", mode: 'link', saveAs: { "${srr}.tsv" }

  input:
  tuple val(srr), path(mqc_fastqc)

  output:
  tuple val(srr), env(status), env(factor)

  """
  # Hypothetical script prints two fields: STATUS and FACTOR
  # STATUS in {low, ok, high}; FACTOR is a float (e.g., 0.5 if high; 1.0 if ok)
  eval_coverage.py $params.genome_size_mb $mqc_fastqc > eval_coverage.tsv
  status=\$(cut -f1 eval_coverage.tsv)
  factor=\$(cut -f2 eval_coverage.tsv)
  echo "Status: \$status // Factor: \$factor"
  """
}


process downsample_reads {
  tag "${srr}_downsample"

  input:
  tuple val(srr), path(r1), path(r2), val(factor)

  output:
  tuple val(srr), path("ds_R1.fq.gz"), path("ds_R2.fq.gz"), emit: reads

  """
  # Check that factor is less than 1.0, otherwise seqtk interprets it differently.
  if (( \$(echo "$factor >= 1" | bc -l) )); then
      echo "[ERROR] Downsample factor ($factor) must be less than 1. Aborting." >&2
      exit 1
  fi

  # sample paired reads deterministically using the same seed for both mates
  seqtk sample -s100 ${r1} ${factor} | gzip -c > ds_R1.fq.gz
  seqtk sample -s100 ${r2} ${factor} | gzip -c > ds_R2.fq.gz
  """
}


process spades {
	tag "${srr}_spades"
	cpus 8
	memory '32 GB'
	publishDir "$params.outdir/assembly", mode: 'link', saveAs: { "${srr}.fna" }

	input:
	tuple val(srr), path(r1), path(r2), path(unpaired), val(use_unpaired)

	output:
	tuple val(srr), path('spades/contigs.fasta'), emit: assembly
	// path 'spades/spades.log' -- if saving this, need to modify publishDir.saveAs

	"""
	if [ "${use_unpaired}" = "true" ]; then SING="-s ${unpaired}"; else SING=""; fi
	spades.py --isolate -t $task.cpus -m 100 --cov-cutoff auto \
		-1 ${r1} -2 ${r2} \$SING -o spades
	"""
}


workflow {
	// main:
	ch_sra_runs = Channel
		.fromPath(params.sra_run_accessions, checkIfExists:true)
		.splitText()
		.map { line -> line.trim() } // trim leading/trailing whitespace
		.filter { it }  // Exclude any empty or whitespace-only lines

	if (params.test) {
		ch_sra_runs = ch_sra_runs.take(2)
	}

	// Prefetch & preprocess
	preproc = ch_sra_runs | prefetch | preprocess_illumina_paired

	// Per-sample MultiQC
	mqc = multiqc(preproc.fastqc) // (srr, mqc_html, mqc_json, mqc_stats)

	// Evaluate coverage using per-sample MultiQC outputs
	decisions = eval_coverage(
	    mqc.map { srr, mqc_html, mqc_json, mqc_fastqc -> tuple(srr, mqc_fastqc) }
	)

	// Join reads with decision by srr
	reads_with_decision = preproc.reads.join(decisions)   // (srr, r1, r2, u, status, factor)

	// Route: too low -> drop; ok -> pass through; high -> downsample (and drop unpaired)
	status_br = reads_with_decision.branch { srr, r1, r2, u, status, factor ->
		ok: status == 'ok'
		high: status == 'high'
		low: true // fallback condition, if status is not ok or high then it should be low
	}

	// Downsample only the 'high' subset
	ds_reads = downsample_reads(status_br.high.map { srr, r1, r2, u, status, factor -> tuple(srr, r1, r2, factor) }).reads

	// Pull together reads for SPAdes (Exclude unpaired for downsampled reads, but we still need to pass them in)
	// (a bit annoying but there's no way to have an "optional input" for spades).
	spades_ok = status_br.ok.map { srr, r1, r2, u, status, factor -> tuple(srr, r1, r2, u, true)}
	unpaired_by_srr = preproc.reads.map { srr, r1, r2, u -> tuple(srr, u) }
	spades_high = ds_reads.join(unpaired_by_srr).map { srr, r1, r2, u -> tuple(srr, r1, r2, u, false)}
	reads_for_spades = spades_ok.mix(spades_high)

	spades(reads_for_spades)
	
	// would be nice to also have a txt of all the decisions 

	/*
	publish:
	sra = sra
	mqc = mqc
	msh = msh
	assemblies = assemblies
	isescan = run_isescan.out
	*/
}

/*
output {
	sra {
		path { srr, sra ->
			"01.sra/${srr}.sra"
		}
	}
	mqc {
		path { srr, report, data_json, fastqc_txt ->
			report >> "02.multiQC/${srr}_report.html"
			data_json >> "02.multiQC/${srr}_data.json"
			fastqc_txt >> "02.multiQC/${srr}_fastqc.txt"
		}
	}
	msh {
		path { srr, msh -> 
			msh >> "04.msh/${srr}.msh"
		}
	}
}
*/