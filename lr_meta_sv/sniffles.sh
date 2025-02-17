#!/bin/bash

#SBATCH --job-name=zurn_snf
#SBATCH --output=zurn_snf_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=asbhatt
#SBATCH --time=8:00:00

#┌────────────┬───────────┬────────────┬───────────────┬───────────────┬───────────────┐
#│ sample     ┆ patientID ┆ collection ┆ reads         ┆ metagenome    ┆ genome        │
#╞════════════╪═══════════╪════════════╪═══════════════╪═══════════════╪═══════════════╡

# To run:
# module load minimap2 seqtk
# ca sniffles (mamba) -- included samtools; if sniffles26 then load it
# sbatch sniffles.sh <sampleID> <reads.fastq.gz> <genome.fasta>

sample="$1"
reads="$2"
genome="$3"

mkdir -p $sample

# Step 0: Index the assembly
chromosome_index="$sample/00.chromosome.mmi"
minimap2 -d $chromosome_index $genome

# Step 1: Align reads to the entire assembly (& process/index sam/bam)
bam_all="$sample/01.all.bam"
minimap2 -t 16 -L -ax map-ont $chromosome_index $reads | \
	samtools view -@ 16 -b | samtools sort -@ 16 > $bam_all

samtools index -@ 16 $bam_all

# Step 2: Filter by MAPQ (e.g., 30, or 10) and contig of interest
bam_filtered="$sample/02.filtered.bam"
samtools view -b -q 10 $bam_all > $bam_filtered

# Step 3: Extract read names and corresponding reads
efm_reads="$sample/03.efm.fastq.gz"
samtools view $bam_filtered | awk '{print $1}' | sort | uniq > $sample/03.efm_reads_list.txt
seqtk subseq $reads $sample/03.efm_reads_list.txt | gzip > $efm_reads

# Step 3.5: Remove intermediate bam files
rm $bam_all
rm $bam_filtered

# Step 4: Re-align & sort the reads
bam_efm="$sample/04.efm.bam"
minimap2 -t 16 -L -ax map-ont $chromosome_index $efm_reads | \
	samtools view -@ 16 -b | samtools sort -@ 16 > $bam_efm

# Note that Sniffles requires sorted and indexed bam files
samtools index -@ 16 $bam_efm

# Step 5: Run Sniffles (v2.2 and v2.6)
snif="$sample/05.sniffles"
sniffles --reference $genome --threads 16 --input $bam_efm --snf $snif.snf --vcf $snif.vcf

# Step 5 alternate: Sniffles v2.6 rerun with mosaic
# snif_mosaic="$sample/05.sniffles.mosaic"
# sniffles --mosaic --reference $genome --threads 16 --input $bam_efm --snf $snif_mosaic.snf --vcf $snif_mosaic.vcf
