#!/usr/bin/env Rscript

# ============================
# IS expansion in Enterococcus analysis
# ============================

# ---------------------- Load packages ----------------------
required_packages <- c(
  "yaml", "ggplot2", "dplyr", "tidyr", "stringr", "reshape2",
  "forcats", "viridis", "ggthemes", "ape", "ggtree", "openxlsx", "tibble", "BiocManager"
)

install_if_missing <- function(pkgs) {
  # Always check for BiocManager first
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing missing package: BiocManager")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing missing package: ", pkg)
      if (pkg == "ggtree") {
		# Install ggtree with BiocManager
        BiocManager::install("ggtree", ask = FALSE, update = FALSE)
      } else {
		# Inatall other packages normally
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
  }
}
install_if_missing(required_packages)

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(reshape2)
  library(forcats)
  library(viridis)
  library(ggthemes)
  library(ape)
  library(ggtree)
  library(openxlsx)
  library(tibble)
})

# ---------------------- Load config ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript is_expansion_phylogeny.R <config.yaml> <distance_threshold>")
}

# ---- Config parameters & args ----
config <- yaml::read_yaml(args[1])
distance_threshold <- as.numeric(args[2])

resolve_path <- function(root, path)
  normalizePath(if (startsWith(path, "/")) path else file.path(root, path))

data_dir <- resolve_path(config$project_root, config$data_dir)
output_dir <- resolve_path(config$project_root, config$output_dir)

# Construct input file paths
gtdbtk_tree_file <- file.path(data_dir, "is_expansion_phylogeny/gtdbtk.bac120.decorated.tree")
dist_matrix_file <- file.path(data_dir, "is_expansion_phylogeny/gtdbtk.bac120.decorated.tree_v2.dist")
lebreton_ref_file <- file.path(data_dir, "is_expansion_phylogeny/representative_genomes.txt")
ref_isescan_file <- file.path(data_dir, "is_expansion_phylogeny/references_isescan_sum.tsv")
isescan_all_file <- file.path(data_dir, "isescan_per_contig.csv")

tree <- read.tree(gtdbtk_tree_file)
dist_matrix <- read.delim(dist_matrix_file, sep = ",", row.names = 1, check.names = F)
ref_genomes <- read.delim(lebreton_ref_file, header = TRUE, stringsAsFactors = FALSE)
ref_isescan <- read.delim(ref_isescan_file, header = TRUE, stringsAsFactors = FALSE)
all_isescan <- read.csv(isescan_all_file, header = TRUE, stringsAsFactors = FALSE)

# prepare output figure paths
fig2A <- file.path(output_dir, "02_A_EFM_ISL3_counts.tsv")
fig2B <- file.path(output_dir, "02_B_IS_family_per_reference.tsv")

# extract leaf labels from tree
tip_nodes <- tree$tip.label

# prepare iTol labels
itol_labels <- data.frame(node_id = tip_nodes)

# Improve formatting of labels for better visualization 
itol_labels <- itol_labels %>%
  dplyr::mutate(label = ifelse(node_id %in% ref_genomes$assembly_accession, node_id, NA)) %>%
  dplyr::mutate(label = ifelse(node_id == "Enterococcus_faecalis_OG1RF_genomic", "Enterococcus faecalis OG1RF", label)) %>%
  dplyr::mutate(label = ifelse(node_id == "Enterococcus_faecalis_V583_genomic", "Enterococcus faecalis V583", label)) %>%
  dplyr::mutate(label = ifelse(node_id == "e_mundtii", "Enterococcus mundtii", label)) %>%
  dplyr::mutate(label = ifelse(node_id == "GCF_000250945.1_ASM25094v1_genomic", "Enterococcus faecium AUS0004", label)) %>%
  dplyr::mutate(label = ifelse(node_id == "AHXU01", "Enterococcus faecium EnGen0043", label)) %>%
  dplyr::mutate(label = ifelse(node_id == "AHWP01", "Enterococcus faecium EnGen0015", label)) %>%
  dplyr::mutate(label = gsub("GCF_[0-9]+\\.[0-9]+_", "", label)) %>% # remove GCF_ and version number
  dplyr::mutate(label = gsub("_genomic", "", label)) %>% # remove _genomic
  dplyr::mutate(label = gsub("_[^_]*$", "", label)) %>% # strip the last underscore and everything after it
  dplyr::mutate(label = gsub("Ente_", "Enterococcus ", label)) %>%
  dplyr::mutate(label = gsub("Tetra_", "Tetragenococcus ", label)) %>%
  dplyr::mutate(label = gsub("Vago_", "Vagococcus ", label)) %>%
  dplyr::mutate(label = gsub("Carnob_", "Carnobacterium ", label)) %>%
  dplyr::mutate(label = gsub("Lacto_", "Lactococcus ", label)) %>%
  dplyr::mutate(label = gsub("Melissococcus_", "Melissococcus ", label)) %>%
  dplyr::mutate(label = gsub("_[^_]*$", "", label)) %>% # strip the last underscore and everything after it
  dplyr::mutate(label = gsub("disp", "dispar", label)) %>%
  dplyr::mutate(label = gsub("cass", "casseliflavus", label)) %>%
  dplyr::mutate(label = gsub("gall", "gallinarum", label)) %>%
  dplyr::mutate(label = gsub("ceco", "cecorum", label)) %>%
  dplyr::mutate(label = gsub("raffi", "raffinosus", label)) %>%
  dplyr::mutate(label = gsub("colu", "columbae", label)) %>%
  dplyr::mutate(label = gsub("pall", "pallens", label)) %>%
  dplyr::mutate(label = gsub("malod", "malodoratus", label)) %>%
  dplyr::mutate(label = gsub("ital", "italicus", label)) %>%
  dplyr::mutate(label = gsub("vill", "villorum", label)) %>%
  dplyr::mutate(label = gsub("asin", "asini", label)) %>%
  dplyr::mutate(label = gsub("phoe", "phoeniculicola", label)) %>%
  dplyr::mutate(label = gsub("sacc", "saccharolyticus", label)) %>%
  dplyr::mutate(label = gsub("haem", "haemoperoxidus", label)) %>%
  dplyr::mutate(label = gsub(" cocc", " caccae", label)) %>%
  dplyr::mutate(label = gsub("gilv", "gilvus", label)) %>%
  dplyr::mutate(label = gsub("mora", "moraviensis", label)) %>%
  dplyr::mutate(label = gsub("sulf", "sulfureus", label)) %>%
  dplyr::mutate(label = gsub("garv", "garvieae", label)) %>%
  dplyr::mutate(label = gsub("malt", "maltaromaticum", label)) %>%
  dplyr::mutate(label = gsub("halo", "halophilus", label)) %>%
  dplyr::mutate(label = ifelse(is.na(label), node_id, label)) %>%
  dplyr::filter(label != "EF013", label != "EF_B_72") # these two have distances from E. faecalis > .15

# using tree distances, get clusters to each representative ---------------

# get references
colnames(ref_genomes) <- "fasta"
representative_genomes <- gsub(".fna", "", ref_genomes$fasta)

# clean distance matrix row and column names, subset distance matrix to references
rownames(dist_matrix) <- gsub(" ", "_", rownames(dist_matrix))
colnames(dist_matrix) <- gsub(" ", "_", colnames(dist_matrix))
dist_matrix <- dist_matrix[ , representative_genomes]

# Convert the distance matrix to long format
dist_long <- dist_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Genome") %>%
  pivot_longer(-Genome, names_to = "Reference", values_to = "Distance")

# filter out EF013 and EF_B_72
dist_long <- dist_long %>%
  dplyr::filter(Genome != "EF013", Genome != "EF_B_72")

# filter out Efm not reaching contiguity thresholds
low_contiguity <- c("FM009", "FM035", "FM_B_05", "FM_B_15")
dist_long <- dist_long %>%
  dplyr::filter(!Genome %in% low_contiguity)

# for each genome, get the reference it is closest to
closest_references <- dist_long %>%
  dplyr::group_by(Genome) %>%
  slice_min(Distance) %>%  # Get the row with the minimum distance
  ungroup()

# Extract the closest reference names
closest_references <- closest_references %>%
  dplyr::filter(Distance <= distance_threshold) %>%
  dplyr::mutate(Assembly = str_extract(Genome, "^[^.]+\\.[0-9]+")) %>%
  dplyr::mutate(Assembly = ifelse(is.na(Assembly), Genome, Assembly)) %>%
  dplyr::mutate(Assembly = gsub("\\.v1", "\\.1", Assembly)) %>%  # Standardize the assembly names
  dplyr::mutate(Assembly = gsub("\\.v2", "\\.2", Assembly)) %>%  # Replace underscores with spaces
  dplyr::left_join(itol_labels, by = c("Reference" = "node_id"))  # Add the class information

# ---- Summarize genome counts ----
total_genomes <- nrow(dist_matrix)
included_genomes <- length(unique(closest_references$Genome))
excluded_genomes <- total_genomes - included_genomes

# 2112, other 16 excluded because of distance filter of 0.01 ^

message(
  sprintf(
    "✅ IS Expansion Analysis: %d genomes included (%d excluded by contiguity or distance threshold ≤ %.3f).",
    included_genomes, excluded_genomes, distance_threshold
  )
)

# get IS counts by family for each genome --------

# prepare IS counts per family for representative genomes
representative_genome_ises_per_family <- ref_isescan %>%
  dplyr::mutate(Sample = gsub(".fna", "", Sample)) %>%
  dplyr::filter(IS_Family != "total") %>%  # Remove the total row
  dplyr::left_join(itol_labels, by = c("Sample" = "node_id")) %>%
  dplyr::group_by(Sample, IS_Family, label) %>%
  dplyr::summarize(count = sum(IS_count), .groups = "drop")

# Define the desired order of label values to match tree
label_order <- c(
  "Lactococcus garvieae", "Carnobacterium maltaromaticum", 
  "Melissococcus plutonius", "Vagococcus lutrae", 
  "Enterococcus faecalis OG1RF", "Enterococcus faecalis V583", 
  "Enterococcus caccae", "Enterococcus moraviensis", 
  "Enterococcus haemoperoxidus", "Enterococcus pallens", 
  "Enterococcus malodoratus", "Enterococcus gilvus", 
  "Enterococcus raffinosus", "Enterococcus avium", 
  "Enterococcus dispar", "Enterococcus asini", 
  "Enterococcus casseliflavus", "Enterococcus gallinarum", 
  "Tetragenococcus halophilus", "Enterococcus saccharolyticus", 
  "Enterococcus italicus", "Enterococcus sulfureus", 
  "Enterococcus cecorum", "Enterococcus columbae", 
  "Enterococcus phoeniculicola", "Enterococcus mundtii", 
  "Enterococcus faecium AUS0004", "Enterococcus faecium EnGen0043", 
  "Enterococcus durans", "Enterococcus hirae", 
  "Enterococcus villorum", "Enterococcus faecium EnGen0015"
)

# Define the desired order of IS_Family values
abundant_is_families <- c("ISL3", "IS30", "IS256", "IS3", "IS6", "IS110")

# Now for all genomes, get IS counts per family ------------
all_ises <- all_isescan %>%
  dplyr::filter(family %in% abundant_is_families) %>%
  dplyr::group_by(sample, family) %>%
  dplyr::summarize(count = sum(count), .groups = "drop")

# filter to the Enterococcus set
# to do so, must convert .1, .2, .3, .4, etc. to _1 and so forth
all_ises_renamed <- all_ises %>%
  dplyr::mutate(sample = gsub("\\.(\\d+)", "_\\1", sample))

# bind representative genome IS data to all IS data
reps_to_bind <- representative_genome_ises_per_family %>%
  dplyr::select(Sample, IS_Family, count) %>%
  dplyr::rename(family = IS_Family, sample = Sample)

all_ises_renamed <- all_ises_renamed %>%
  dplyr::bind_rows(reps_to_bind)

all_ises_per_family <- closest_references %>%
  dplyr::left_join(all_ises_renamed, by = c("Genome" = "sample"))

# 1) list of genomes for each label
genomes_per_label <- all_ises_per_family %>%
  distinct(label, Genome)

# 2) list of families you want to include
families <- unique(all_ises_per_family$family)

# 3) full grid: every (label, Genome, family)
full_grid <- genomes_per_label %>%
  tidyr::crossing(family = families)

# 4) bring in observed counts, fill missing with 0
full_counts <- full_grid %>%
  dplyr::left_join(
    all_ises_per_family %>% dplyr::select(label, Genome, family, count),
    by = c("label", "Genome", "family")
  ) %>%
  dplyr::mutate(count = tidyr::replace_na(count, 0))

# 5) summarize per label x family
plot_data <- full_counts %>%
  dplyr::filter(family %in% abundant_is_families) %>%
  group_by(label, family) %>%
  dplyr::summarise(
    average = mean(count),
    median = median(count),
    stdev = as.numeric(if (n() > 1) sd(count) else 0),    # sd on single sample -> 0
    max_count = max(count),
    # if all zeros, set max_genome to NA; otherwise pick genome with max count
    max_genome = if (all(count == 0)) NA_character_ else Genome[which.max(count)][1],
    n = n(),   # number of genomes per label (should match genome count)
    .groups = "drop"
  ) %>%
  left_join(itol_labels, by = "label") %>%
  dplyr::mutate(
    label = factor(label, levels = label_order),
    family = factor(family, levels = abundant_is_families)
  ) %>%
  dplyr::mutate(average = as.numeric(sprintf("%.3f", average)))

# pivot plot_data
all_genome_ises_per_family_wide <- plot_data %>%
  dplyr::select(node_id, family, average) %>%
  dplyr::mutate(average = as.numeric(average)) %>%
  dplyr::filter(!is.na(node_id)) %>%
  pivot_wider(names_from = family, 
              values_from = average)

all_genome_ises_abundant_family_wide <- plot_data %>%
  dplyr::select(node_id, family, average) %>%
  dplyr::mutate(average = as.numeric(average)) %>%
  dplyr::filter(!is.na(node_id)) %>%
  pivot_wider(names_from = family, 
              values_from = average) %>%
  dplyr::mutate(other = rowSums(dplyr::select(., -node_id, -all_of(abundant_is_families[-length(abundant_is_families)])), na.rm = TRUE)) %>%
  dplyr::select(node_id, any_of(abundant_is_families))

# adjust labels for readability
fig2B_output <- all_genome_ises_abundant_family_wide %>%
  rename(label = node_id) %>%
  mutate(
    # ---- Clean up GenBank junk and underscores ----
    label = str_replace(label, "^GCF_\\d+\\.\\d+_", ""),          # remove GCF prefix
    label = str_replace(label, "_ASM\\d+v\\d+", ""),              # remove assembly ID
    label = str_replace(label, "_genomic$", ""),                   # remove suffix
    label = str_replace_all(label, "_", " "),                      # underscores → spaces
    label = str_trim(label),
    
    # ---- Expand shorthand labels ----
    label = case_when(
      label == "e mundtii" ~ "Enterococcus mundtii",
      label == "AHWP01" ~ "Enterococcus lactis (B)",
      label == "AHXU01" ~ "Enterococcus faecium (A2)",
      label == "ASM25094v1" ~ "Enterococcus faecium (A1)",
      TRUE ~ label
    ),
    
    # ---- Standardize genus names ----
    label = str_replace(label, "^Ente ", "Enterococcus "),
    label = str_replace(label, "^Lacto ", "Lactococcus "),
    label = str_replace(label, "^Carnob ", "Carnobacterium "),
    label = str_replace(label, "^Tetra ", "Tetragenococcus "),
    label = str_replace(label, "^Vago ", "Vagococcus "),
    label = str_replace(label, "^Melissococcus ", "Melissococcus "),
    
    # ---- Abbreviate genus to initial + period ----
    label = str_replace(label, "^Enterococcus ", "E. "),
    label = str_replace(label, "^Lactococcus ", "L. "),
    label = str_replace(label, "^Carnobacterium ", "C. "),
    label = str_replace(label, "^Tetragenococcus ", "T. "),
    label = str_replace(label, "^Vagococcus ", "V. "),
    label = str_replace(label, "^Melissococcus ", "M. "),
    
    # ---- Custom display names ----
    label = case_when(
      label == "E. faecalis V583" ~ "E. faecalis (clinical)",
      label == "E. faecalis OG1RF" ~ "E. faecalis (lab)",
      label == "E. faecium EnGen0015" ~ "E. lactis (B)",
      TRUE ~ label
    )
  )

# Save the wide format data
write.table(fig2B_output, fig2B, sep = "\t", quote = FALSE, row.names = FALSE)
message("✅ IS family counts per reference saved to: ", fig2B)

# get genomes for clade B, clade A1, clade A2 ------------------------
efm_ranges <- closest_references %>%
  dplyr::filter(grepl("faecium", label)) %>%
  dplyr::mutate(Clade = case_when(
    grepl("aus0004", label, ignore.case = T) ~ "clade A1",
    grepl("engen0043", label, ignore.case = T) ~ "clade A2",
    grepl("engen0015", label, ignore.case = T) ~ "clade B",
    TRUE ~ "other"
  )) %>%
  dplyr::select(Genome, Clade)

# print counts per clade --------------------------------------------------
clade_counts <- efm_ranges %>%
  dplyr::count(Clade, name = "n") %>%
  dplyr::arrange(desc(n))

message("📊 E. faecium genomes by clade:")
print(clade_counts)

# write ISL3 counts for Efm genomes to file -------------------------------
isl3_counts <- efm_ranges %>%
  # start with all Efm genomes
  dplyr::select(Genome, Clade) %>%
  # join the ISL3 counts (may have missing genomes)
  dplyr::left_join(
    all_ises_per_family %>%
      dplyr::filter(family == "ISL3") %>%
      dplyr::select(Genome, count),
    by = "Genome"
  ) %>%
  # replace NA counts with 0
  dplyr::mutate(count = tidyr::replace_na(count, 0))

write.table(isl3_counts, fig2A, sep = "\t", quote = FALSE, row.names = FALSE)
message("✅ ISL3 counts for Efm genomes saved to: ", fig2A)
