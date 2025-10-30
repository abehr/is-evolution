import 'common.just'

default:
  @echo "ISL3 proliferation in Enterococcus faecium"
  @just --list

#------------------------------------------------------------------------------
# 🧱 Setup recipes
#------------------------------------------------------------------------------

install:
  @python3 -c "import sys; assert sys.version_info >= (3,10), '❌ Python >=3.10 required, found: ' + sys.version.split()[0]"
  @python3 -c "import sys; v = sys.version_info; import sys; sys.stdout.write('⚠️ Python version is below 3.13 (tested version); this should work but has not been tested as extensively.\n' if (v >= (3,10) and v < (3,13)) else '')"
  @echo "✅ Python version OK: $$(python3 --version | cut -d' ' -f2)"
  @echo "Creating Python virtual environment..."
  python3 -m venv .venv
  source .venv/bin/activate && pip install -U pip wheel setuptools
  source .venv/bin/activate && pip install -e .[all]
  cp -n config.example.yaml ../config.yaml || true
  @echo "✅ Environment ready. Edit ../config.yaml before running."

install-bare:
  @echo "Setting up environment..."
  {{python}} -m pip install -U -e .
  cp -n config.example.yaml ../config.yaml || true
  @echo "✅ Environment ready. Edit ../config.yaml before running."

cleanup:
  @echo "🧹 Cleaning up Python virtual environment and build artifacts..."
  # pip-uninstall in case we're outside of a venv
  @echo "↪ Uninstalling package(s) from current Python (if present)..."
  {{python}} -m pip uninstall -y is-evolution || true
  rm -rf .venv
  find . -type d -name "__pycache__" -exec rm -rf {} +
  find . -type d -name ".pytest_cache" -exec rm -rf {} +
  find . -type d -name "*.egg-info" -exec rm -rf {} +
  @echo "✅ Cleanup complete. Run just setup-venv to rebuild."


#------------------------------------------------------------------------------
# 📦 Example data processing steps
#------------------------------------------------------------------------------

# Requires ncbi-datasets-cli (can be conda-installed)
download-ncbi-complete taxon:
  # Optionally, could install it here:
  # micromamba install -c conda-forge ncbi-datasets-cli
  @echo "Download {{taxon}} complete genomes"
  which datasets # make sure datasets is installed
  datasets download genome taxon {{taxon}} \
    --assembly-level complete \
    --include genome,gff3 \
    --dehydrated \
    --filename {{taxon}}.zip


# build-datasets does not actually need to be run, as the outputs of this have already been produced
# and can be downloaded. However, we provide these helper scripts and this example routine
# so it's easy to see how they were generated, and because they could also be modified to 
# run on other taxa/datasets, though the corresponding intermediate files (such as concatenated 
# geNomad and ISEScan results) would need to be produced as well. Additional configuration and
# files would be needed in order to directly run this recipe. 
build-datasets:
  # Parse NCBI complete genomes info from downloaded NCBI datasets results
  # ↘ source.ncbi_complete_datasets
  # ↗ data.ncbi.complete_genomes
  {{python}} dataproc/ncbi_complete.py

  # Combine NCBI complete genomes & NCBI long-read assemblies, and apply most QC-filters.
  # ↘ data.ncbi.complete_genomes, data.ncbi.nanopore_samples
  # ↗ data.ncbi.combined_qc_filtered
  {{python}} dataproc/combine_all_ncbi.py

  # Filter to samples with paired geNomad and ISEScan results, and only spp with >10 samples
  # ↘ data.ncbi.combined_qc_filtered
  # ↗ data.ncbi.combined_final_samples
  {{python}} dataproc/final_samples.py

table-is-counts:
  @echo "[Supplementary Data 2] NCBI-complete & NCBI-longread sample info (genome length, contig info, IS counts)"
  @echo "[Supplementary Data 3] SHC-iso sample info (genome length, IS counts where necessary for analysis)"
  {{python}} dataproc/table_is_counts.py


#--------------------------------------------------------------------------------------
# 🔍 Steps for curation of short-read datasets for IS historical timelines
#--------------------------------------------------------------------------------------

# Note that the final output(s) of these recipes are available in our Zenodo project,
# so these steps do not need to be run, but we provide them here to show how the 
# short read dataset accessions were curated. 

# Query NCBI API to get info on all relevant datasets & metadata for a given organism
entrez-esearch-fetch wdir ncbi_auth_email ncbi_api_key organism:
  {{python}} {{repo_root}}/src/biobehr/ncbi_entrez_direct.py \
    --email {{ncbi_auth_email}} --key {{ncbi_api_key}} \
    -o {{wdir}} \
    --organism '{{organism}}' \
    --platform illumina \
    --source genomic \
    --strategy wgs \

# Some semi-generalized preprocessing of results (to prep for more ad-hoc/specific filtering)
entrez-results-preprocess wdir genome_mb:
  {{python}} process_entrez_results.py {{genome_mb}} \
    {{wdir}}/samples.jsonl \
    {{wdir}}/experiments.jsonl \
    {{wdir}}/processed_entrez_results.csv


### E. faecium short-read dataset curation
efaecium-get-shortread-accessions wdir ncbi_auth_email ncbi_api_key:
  just analysis/03_shortread_is_estimate_timeline/entrez-esearch-fetch {{wdir}} {{ncbi_auth_email}} {{ncbi_api_key}} 'Enterococcus faecium'
  just analysis/03_shortread_is_estimate_timeline/entrez-results-preprocess {{wdir}} 2.9
  {{python}} efm_filter_sra.py {{wdir}}/processed_entrez_results.csv {{wdir}}/efm_filtered_runs.txt


### S. epidermidis short-read datset curation
staphepi-get-shortread-accessions wdir ncbi_auth_email ncbi_api_key:
  just analysis/03_shortread_is_estimate_timeline/entrez-esearch-fetch {{wdir}} {{ncbi_auth_email}} {{ncbi_api_key}} 'Staphylococcus epidermidis'
  just analysis/03_shortread_is_estimate_timeline/entrez-results-preprocess {{wdir}} 2.5
  {{python}} analysis/03_shortread_is_estimate_timeline/epi_filter_sra.py {{wdir}}/processed_entrez_results.csv {{wdir}}/epi_filtered_runs.txt


# NOTE for above: wdir must be provided as an absolute path since we're calling sub-justfiles and passing it in.


#------------------------------------------------------------------------------
# 📊 Analyses & generating output
#------------------------------------------------------------------------------

# Generate plots for Figure 1 and related Extended Data Figs
pathogen-is-survey:
  @echo "[Fig. 1B] IS counts per ESKAPEE organism"
  {{python}} analysis/01_pathogen_is_survey/plot_eskape_is_counts.py
  
  @echo -e "\n[Extended Data Fig. 1] IS counts/spp across ESKAPEE genera"
  {{python}} analysis/01_pathogen_is_survey/plot_is_per_taxon.py

  @echo -e "\n[Fig. 1C][Extended Data Fig. 2] IS counts/spp (distribution across chromosome vs. plasmids)"
  {{python}} analysis/01_pathogen_is_survey/plot_taxon_is_per_contig_type.py


# Example of how the pairwise ANI cluster analysis can be invoked
pairwise-ani-clusters:
  # NOTE: requires additional configuration in config.yaml:
  # source.mash_pairwise_ani fields must be set; these files are generated from public data.
  # Mash pairwise distance matrix from the publicly-available genome assemblies,
  # as specified in the manuscript Methods and Extended Data Fig. 3 legend. 
  # The distance matrix is not included in our data upload due to its large 
  # file size and because it can be recreated from public data.
  @echo "[Extended Data Fig. 3] Genomic diversity per ESKAPEE organism and IS copy range per cluster"
  {{r}} analysis/01_pathogen_is_survey/pairwise_ani.R $PROJECT_CONFIG 0.01 25


# Plot IS tpase sequence conservation/diversity & define flank lengths for common IS families.
tpase-seqs-flanks-boundaries:
  # NOTE: requires downloading additional public data (E. faecium genomes & GFFs from NCBI) to run. 
  # NOTE: requires MAFFT (or Muscle5 with a one-line change).

  @echo "[Fig. 1D] Plot IS tpase sequence conservation/diversity; output common tpase seqs & sequence context."
  {{python}} analysis/01_seqs_flanks_boundaries/exact_copy_analysis.py

  @echo -e "\n[Data for Extended Data Fig. 5] Analyze context of common tpase in order to determine entropy-based IS boundaries."
  {{python}} analysis/01_seqs_flanks_boundaries/calculate_is_flank_boundaries.py


# Plots for Figure 2 and related extended data figs
isl3-taxonomic-distribution:
  @echo "[Data for Fig. 2A-B] ISL3 expansion phylogeny"
  {{r}} analysis/02_isl3_taxonomic_distribution/is_expansion_phylogeny.R $PROJECT_CONFIG 0.01

  @echo -e "\n[Fig. 2C] ISL3 counts/spp across gram-positive cocci"
  {{python}} analysis/02_isl3_taxonomic_distribution/plot_isl3_per_cocci.py


### Main plots for Figure 3

# `name` refers to the name of a directory within `shortread_is_timeline` containing 
# concatenated results after running the pre.nf and post.nf pipelines and performing QC.
# Included in the Zenodo are "efaecium" and "staph_epi"
efm-is-timeline:
  echo "Estimate E. faecium IS counts from ISEScan and coverage info in short-read datasets"
  {{python}} analysis/03_shortread_is_estimate_timeline/count_ise.py efaecium

  echo -e "\nPlot E. faecium IS counts over time"
  {{python}} analysis/03_shortread_is_estimate_timeline/efm_plot_ise.py


epi-is-timeline:
  echo "Estimate S. epidermidis IS counts from ISEScan and coverage info in short-read datasets"
  {{python}} analysis/03_shortread_is_estimate_timeline/count_ise.py staph_epi
  
  echo -e "\nPlot S. epidermidis IS counts over time"
  {{python}} analysis/03_shortread_is_estimate_timeline/epi_plot_ise.py



# Example of how the Pangraph analysis workflow can be invoked; arguments should be modified based on your workflow.
pangraph-structural-variation dataset='':
  @echo "[Data for Fig. 4C-D] Pangraph junction statistics (IS-mediated structural variation)"
  @if [ -z "{{dataset}}" ]; then \
    {{python}} analysis/04_pangraph_structural_variation/junction_stats.py -h; \
  else \
    {{python}} analysis/04_pangraph_structural_variation/junction_stats.py \
      ../local/pangraph/{{dataset}}/pangraph_output \
      ../local/pangraph/{{dataset}}/pangraph_input_fasta \
      ../local/pangraph/{{dataset}}/gff \
      ../output/pangraph/{{dataset}}; \
  fi

longitudinal-structural-variation:
  @echo "[Fig. 5B] Within-sample structural variation"
  {{python}} analysis/05_longitudinal_metagenomic/sniffles_sv_within_sample.py


rnaseq-diff-expression:
  # The RNA-Seq differential expression analysis takes in the following params:
  # (1) Project config.yaml
  # (2) Feature count matrix file. These data are included in our manuscript 
  #     in Supplementary Data 7, but must be converted from Excel to plaintext
  #     format to run this R script.
  # (3) Reference GFF file. This can be generated by running Bakta on the 
  #     p11568 t4 assembled chromosome.
  # (4) Significance threshold (1e-25 was used in our manuscript).

  # NOTE: Certain installations of R (especially macOS via Homebrew) may run into 
  # errors while attempting to install this script's dependencies, such as the DESEq2
  # library, which has some system-level dependencies. If this occurs, it might be
  # best to switch to a different R environment.
  @echo "[Data for Fig. 6B-C][Extended Data Fig. 11] RNA-Seq differential expression analysis \
  (ISL3 variant effects on neighboring gene expression)"
  {{r}} analysis/06_transcriptomic/rnaseq_analysis.R $PROJECT_CONFIG \
    ../local/rna_seq/multicov_counts_v2.txt \
    ../local/rna_seq/p11568_04.gff \
    1e-25


genes-adjacent-is:
  # NOTE: requires additional configuration in config.yaml:
  # - source.poppunk_tree_genomes (list of desired sample IDs to analyze)
  # - source.gff (directory containing GFFs to analyze).
  # GFFs can be downloaded from NCBI (for public data) or generated
  # by running Bakta (for SHC data on Zenodo).
  @echo "[Supplementary Data 8] ISL3->FolT and other genes adjacent to IS transposases in E. faecium"
  {{python}} analysis/06_isl3_proximal_genes/find_is_proximal_genes.py



