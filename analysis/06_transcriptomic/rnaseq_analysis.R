#!/usr/bin/env Rscript

# ============================
# RNA-seq analysis pipeline, analysis related to Figure 6A-C 
# and Extended Data Fig 11
# ============================

# ---------------------- Load config & parameters ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript rnaseq_analysis_publish.R <config.yaml> <counts_file> <gff_file> <significance_threshold>")
}

config <- yaml::read_yaml(args[1])

output_dir <- normalizePath(if (startsWith(config$output_dir, "/")) config$output_dir else file.path(config$project_root, config$output_dir))

counts_file <- args[2]
gff_file <- args[3]
significance_threshold <- as.numeric(args[4])

# ---- Package setup ----
required_packages <- c(
  "yaml", "ggplot2", "dplyr", "tidyr", "viridis", "ggrepel",
  "DESeq2", "EnhancedVolcano", "ape", "stringr", "optparse"
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
      if (pkg == "DESeq2") {
        BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
      } else if (pkg == "EnhancedVolcano") {
		BiocManager::install("EnhancedVolcano", ask = FALSE, update = FALSE)
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
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(ggrepel)
  library(DESeq2)
  library(EnhancedVolcano)
  library(ape)
  library(stringr)
})

# ---- Directory setup ----
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("Reading input data...")

# ---- Load input data ----
df <- read.table(counts_file, sep = "\t", header = FALSE, fill = TRUE, quote = "")
colnames(df) <- c(
  "seqid", "source", "type", "start", "end", "score", "strand", "phase",
  "attributes", "t1_1", "t1_2", "t1_3", "t1_4", "t4_1", "t4_2", "t4_3", "t4_4"
)

counts <- df %>%
  separate(attributes, sep = ";", into = c("ID","Name","locus_tag","product","dbxref"),
           extra = "drop", fill = "right") %>%
  dplyr::select(-dbxref) %>%
  dplyr::mutate(product = str_remove(product, 'product=')) %>%
  dplyr::mutate(Name = str_remove(Name, 'Name=')) %>%
  dplyr::mutate(locus_tag = str_remove(locus_tag, 'locus_tag=')) %>%
  dplyr::mutate(ID = str_remove(ID, 'ID=')) %>%
  separate(seqid, sep = "_NODE", into = c('sample','rest'), remove = FALSE,
           extra = "drop", fill = "right") %>%
  dplyr::select(-rest, -sample)

# ---- Load GFF annotations ----
features <- ape::read.gff(gff_file)
annotations <- features %>%
  dplyr::filter(type != "region") %>%
  dplyr::filter(seqid %in% c("contig_1", "contig_2", "contig_4")) %>%
  separate(attributes, sep = ";", into = c("ID","Name","locus_tag","product","dbxref"),
           extra = "drop", fill = "right") %>%
  separate(dbxref, sep = "UniRef:", into = c("first","UniRef_1","UniRef_2","UniRef_3"), 
           extra = "drop", fill = "right") %>%
  separate(first, sep = ',', into = c('att1','att2','att3','att4','att5'), 
           extra = "drop", fill = "right") %>%
  dplyr::mutate(UniRef_1 = str_remove(UniRef_1, ',')) %>%
  dplyr::mutate(UniRef_2 = str_remove(UniRef_2, ',')) %>%
  dplyr::mutate(product = str_remove(product, 'product=')) %>%
  dplyr::mutate(Name = str_remove(Name, 'Name=')) %>%
  dplyr::mutate(locus_tag = str_remove(locus_tag, 'locus_tag=')) %>%
  dplyr::mutate(ID = str_remove(ID, 'ID=')) %>%
  dplyr::select(-score) 

# ==========================================================
# INTERNAL PARAMETERS
# ==========================================================
sample_names <- c("t1_1", "t1_2", "t1_3", "t4_1", "t4_2", "t4_3", "t4_4")
conditions <- c("wt", "wt", "wt", "is", "is", "is", "is")
reference_condition <- "wt"
test_condition <- "is"
folT_locus <- c("OAJOIO_14355", "OAJOIO_14350", "OAJOIO_14345")
agr_operon <- c("OAJOIO_01100", "OAJOIO_01095", "OAJOIO_01090", "OAJOIO_01085", "OAJOIO_01080", "OAJOIO_01075")

# ==========================================================
# DESEQ2 ANALYSIS
# ==========================================================

# chromosomal contigs only - as determined by genomad
deseq <- counts %>%
  filter(type != "region", seqid %in% c("contig_1", "contig_2", "contig_4")) %>%
  select(ID, starts_with("t")) %>%
  select(-type, -t1_4) # remove outlier due to poor alignment fraction

tags <- deseq$ID
count_data <- deseq %>%
  select(-ID) %>%
  mutate_all(as.numeric) %>%
  as.matrix()
rownames(count_data) <- tags

colData <- data.frame(row.names = sample_names, condition = factor(conditions))
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~condition)
dds$condition <- relevel(dds$condition, ref = reference_condition)
dds <- estimateSizeFactors(dds)

# Initial number of genes
n_genes_initial <- nrow(dds)

# remove genes that are strictly 0 counts in either condition
zero_counts <- rowSums(counts(dds, normalized = TRUE)[, colData$condition == reference_condition]) == 0 |
  rowSums(counts(dds, normalized = TRUE)[, colData$condition == test_condition]) == 0
dds <- dds[!zero_counts, ]
n_genes_after_zero <- nrow(dds)

# keep genes with at least 5 counts in at least 3 samples
dds <- dds[rowSums(counts(dds, normalized = TRUE) >= 5) > 3, ]
n_genes_after_low <- nrow(dds)

cat("Number of genes at each filtering step:\n")
cat("  Before filtering:", n_genes_initial, "\n")
cat("  After removing zero-count genes:", n_genes_after_zero, "\n")
cat("  After removing low-count genes:", n_genes_after_low, "\n\n")

dds <- DESeq(dds)
res <- results(dds)
res$Name <- rownames(res)

# ---- Merge results ----
res_df <- as.data.frame(res) %>%
  left_join(counts, by = c("Name" = "ID"))
write.csv(res_df, file.path(output_dir, "DESeqResults.csv"), quote = FALSE, row.names = FALSE)

# ---- Identify DE genes ----
de_genes <- res[res$padj < significance_threshold & !is.na(res$padj), ]
upregulated <- sum(de_genes$log2FoldChange > 2)
downregulated <- sum(de_genes$log2FoldChange < -2)
cat("Summary of DE Genes:\n",
    "Total:", nrow(de_genes), "\n",
    "Upregulated:", upregulated, "\n",
    "Downregulated:", downregulated, "\n")

# ==========================================================
# PLOTTING
# ==========================================================
keyvals <- case_when(
  res$Name %in% folT_locus ~ 'green',
  res$Name %in% agr_operon ~ 'gold',
  abs(res$log2FoldChange) > 2 & res$padj < significance_threshold ~ 'red',
  TRUE ~ 'gray'
)
names(keyvals) <- case_when(
  keyvals == 'red' ~ 'DE',
  keyvals == 'gold' ~ 'agr operon',
  keyvals == 'green' ~ 'folT',
  TRUE ~ 'Other'
)

p <- EnhancedVolcano(
  res,
  lab = NA,
  x = 'log2FoldChange',
  y = 'pvalue',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'p-value'),
  colCustom = keyvals,
  pCutoff = significance_threshold,
  FCcutoff = 2.0,
  pointSize = 2.0,
  title = 'RNA-seq Differential Expression',
  subtitle = 'All features',
  legendPosition = 'right'
)

ggsave(file.path(output_dir, "Extended_Data_Fig_11A_Volcano.pdf"), p, width = 8, height = 6)
message("✅ Analysis complete. Results saved to: ", output_dir)