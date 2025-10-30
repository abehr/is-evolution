# ISL3 elements in _Enterococcus faecium_

This repository contains the code, data, and workflows underlying the analyses presented in [our paper](#Citation). 

We examine:
- the landscape of insertion sequences (IS elements) in modern ESKAPEE pathogens,
- the recent proliferation of ISL3 in _E. faecium_,
- extensive structural variation in hospital-adapted lineages via ISL3 activity,
- IS-mediated structural variation in the gut microbiome of hospital patients via long-read metagenomics, and
- ISL3-driven evolution in _E. faecium_ and its potential role in pathogenic adaptation.

> [!WARNING]
> This work is currently under review; code and manuscript may change.

## Quick start

```bash
# Check that required tools are installed
python --version; Rscript --version; just --version

# Create and enter a new project folder
mkdir isl3_efaecium_analysis && cd $_

# Clone this repository
git clone https://github.com/abehr/is-evolution.git

# Unpack preprocessed data (after downloading from Zenodo)
unzip analyses_preprocessed_source_data.zip -d data

# Install package and dependencies
cd is-evolution
just install

# Configure project_root, data_dir, and output_dir
vi ../config.yaml

# Example: generate IS count tables (Supp. Data 2-3)
just table-is-counts

# Example: plot ESKAPEE IS survey results (Fig. 1 + Extended)
just pathogen-is-survey
```


## Repository structure

- **dataproc/** - Scripts to preprocess raw data.
  
  Their outputs are archived on Zenodo and used in downstream analyses.

- **analysis/** - Scripts that generate the figures and tables for the study.
  
  Subfolders are prefixed by the corresponding figure or section number.

- **workflows/** - Nextflow pipelines and configuration examples for producing the input data used in analyses.
  
  These were tuned for our HPC environment and may require adaptation for your system.

- **src/** - Shared code modules.
	+ `src/biobehr`: general utilities (e.g. GFF parsing, NCBI API access)
	+ `src/efm`: project-specific helpers (e.g. configuration handling, plot styling)

## Justfile and analysis recipes

The root-level `justfile` defines reproducible analysis recipes.
Each recipe documents how we generated a given figure or dataset's source material.

For example:

- `pathogen-is-survey` runs all analyses in `analysis/01_pathogen_survey`, producing Figs. 1B–C and Extended Data Figs. 1–2.

	This can be run directly using the preprocessed Zenodo data.
- `pangraph-structural-variation` demonstrates how the PanGraph analysis was configured for Fig. 4 and Supplementary Data 4.

	This requires additional local data and configuration, but you can run the recipe with no arguments to see a more detailed explanation.

While the `just` recipes provide more standardized and reproducible commands, the underlying scripts can also be run directly. 


## Installation

### Prerequisites

- Python ≥ 3.10 (3.13 recommended)
- R ≥ 4.4
- [just](https://github.com/casey/just) (command runner)

We recommend installing via a package manager. Tested environments include macOS (Homebrew) and macOS/Linux (Micromamba).

```bash
# Homebrew
brew install python3 just r

# Conda (existing env)
micromamba install -c conda-forge python=3.13 just

# Conda (new env)
micromamba create -n is_evo_env -c conda-forge python=3.13 just
micromamba activate is_evo_env
```

> **Note:** Google Chrome is required for SVG export in plotting scripts. 
> These scripts are lightweight and optimized for laptop use, but you can comment out the SVG-export lines if running in a non-GUI environment (e.g. over SSH).


### Data download and setup

Preprocessed source data (outputs from `dataproc/` and concatenated summary tables) can be downloaded from our Zenodo archive.

_(Public link to be added upon publication.)_

These files enable analyses to be run locally without reprocessing raw data.

> **Note:** This does not include genomic data that are already publicly available (such as genomes or GFF annotation files that are available on NCBI/SRA). Therefore, while some of our analyses can be run "out of the box", others require additional data download/processing and configuration.

Example project layout:

```bash
efm_isl3_analysis/  # project root
├── is-evolution/   # cloned repository
└── data/           # Zenodo data (CSV and other files)
```

To verify setup and initialize the environment:
```bash
cd is-evolution
just                # list available recipes
just install        # setup project & dependencies using a venv; create config.yaml
```

Alternatively, if you prefer to manage dependencies yourself (e.g. via a conda env), you can run `just install-bare` to only setup the package (without a virtual env & dependencies) and copy the example config.

After either installation method, edit the new `../config.yaml` to set `project_root`, `data_dir`, and `output_dir` before running analyses.
Relative paths are resolved from `project_root`; absolute paths can also be used.

> **Note:** Config location can be updated via the `.env` file (automatically loaded by the `justfile`).

## Running analyses

Many analyses can be launched directly via `just`, for example:

```bash
just pathogen-is-survey
```

This command executes all of the main scripts in `analysis/01_pathogen_survey` producing Figs. 1B-C and Extended Data Fig. 1-2.

Some other recipes require additional data (from our Zenodo archive, our NCBI project, or public NCBI sources) or preprocessing steps.


## Analysis map

This table shows where each of the main analyses lives, which scripts comprise it, and how it relates to figures/tables. For runnable commands, see the corresponding recipe in the justfile. 

| Folder | Main scripts | Outputs | `just` recipe |
|---|---|-----|------|
| **01_pathogen_is_survey/** | `plot_eskape_is_counts.py`, `plot_is_per_taxon.py`, `plot_taxon_is_per_contig_type.py` | Fig. 1B–C, Ext. Data Figs. 1–2 | `pathogen-is-survey` |
| **01_seqs_flanks_boundaries/** | `exact_copy_analysis.py`, `calculate_is_flank_boundaries.py` | Fig. 1D, data for Ext. Data Fig. 5 | `tpase-seqs-flanks-boundaries` |
| **02_isl3_taxonomic_distribution/** | `is_expansion_phylogeny.R`, `plot_isl3_per_cocci.py` | Data for Fig. 2A–B, Fig. 2C | `isl3-taxonomic-distribution` |
| **03_shortread_is_estimate_timeline/** | `count_ise.py`, `ise_plots.py` | Fig. 3C, Ext. Data Fig. 6C | `efm-is-timeline` |
| **04_pangraph_structural_variation/** | `junction_stats.py` | Data for Fig. 4C–D | `pangraph-structural-variation` |
| **05_longitudinal_metagenomic/** | `sniffles_sv_within_sample.py` | Fig. 5B | `longitudinal-structural-variation` |
| **06_transcriptomic/** | `rnaseq_analysis.R` | Fig. 6B–C, Ext. Data Fig 11 | `rnaseq-diff-expression` |
| **06_isl3_proximal_genes/** | `find_is_proximal_genes.py` | Supp. Data 8 | `genes-adjacent-is` |



## Citation

If you use any of the code found in this repository, please cite our [preprint](https://www.biorxiv.org/content/10.1101/2025.03.16.643550v1).

**Replicative selfish genetic elements are driving rapid pathogenic adaptation of _Enterococcus faecium_** \

Matthew Grieshop & Aaron Behr et al., _bioRxiv_ (2025) \
doi: https://doi.org/10.1101/2025.03.16.643550
