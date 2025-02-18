#!/bin/bash

#SBATCH --job-name=zurn_snf
#SBATCH --output=zurn_snf_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --account=asbhatt
#SBATCH --time=12:00:00

#┌────────────┬───────────┬────────────┬───────────────┬───────────────┬───────────────┐
#│ sample     ┆ patientID ┆ collection ┆ reads         ┆ metagenome    ┆ genome        │
#╞════════════╪═══════════╪════════════╪═══════════════╪═══════════════╪═══════════════╡

# To run: ca sniffles_ngmlr (mamba)
# Sniffles v2.6.0
# ngmlr v0.2.7

# patient="$1"
sample="$2" # Current sample+reference name
ref="$3" # Reference sample contigs fasta (path relative to patient dir)
reads="$4" # Raw reads for current sample (path relative to patient dir)

if [ ! -f "$ref" ] || [ ! -f "$reads" ]; then
	echo "Error: One or more required files/directories are missing."
	exit 1
fi

if [ -f "$sample.bam" ]; then
	echo "Alignment file already exists; no need to do anything."
	exit 0
fi

module load samtools
if ! command -v sniffles &> /dev/null || ! command -v ngmlr &> /dev/null || ! command -v samtools &> /dev/null; then
	echo "Error: Required command(s) missing."
	exit 1
fi

# Run NGMLR
ngmlr -r $ref -q $reads -x ont -t 16 --bam-fix -o $sample.sam

# Fix sam format
awk 'BEGIN { OFS="\t" } 
     /^@/ { print; next } 
     { if ($5 < 0) $5 = 60; print }' $sample.sam > $sample.fixed.sam


# Index, sort, etc.
samtools view -b $sample.fixed.sam > $sample.bam
samtools sort -@ 16 -o $sample.sorted.bam $sample.bam
samtools index -@ 16 $sample.sorted.bam
rm $sample.sam
rm $sample.fixed.sam
rm $sample.bam

# samtools view -@ 16 -b $sample.fixed.sam | samtools sort -@ 16 | samtools index -@ 16

# Run Sniffles v2.6, normal and mosaic
snif="$sample.sniffles"
sniffles --reference $ref --threads 16 --input $sample.sorted.bam --minsvlen 100 --max-splits-base 10 --max-splits-kb 0.5 --vcf $snif.vcf
sniffles --reference $ref --threads 16 --input $sample.sorted.bam --minsvlen 100 --max-splits-base 10 --max-splits-kb 0.5 --vcf $snif.mosaic.vcf --mosaic 
