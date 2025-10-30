#!/usr/bin/env Rscript

# ---------------------- Libraries ----------------------
suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(igraph)
})

# ---------------------- Load config ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript pairwise_ani.R <config.yaml> <distance_threshold> <cluster_min_size>")
}

config <- yaml::read_yaml(args[1])
threshold <- as.numeric(args[2]) # distance threshold
minsize <- as.numeric(args[3]) # cluster min size

# ---------------------- Input data ----------------------
mash_file <- config$source$mash_pairwise_ani$mash_dists_file
distance_dir <- config$source$mash_pairwise_ani$mash_dist_dir

# --- Input existence checks ---
if (!file.exists(mash_file) || dir.exists(mash_file)) {
  stop("Error: mash_file does not exist or is not a regular file: ", mash_file)
}
if (!dir.exists(distance_dir) || file.exists(distance_dir)) {
  stop("Error: distance_dir does not exist or is not a directory: ", distance_dir)
}

# ---------------------- Output files --------------------
density_plot_fp <- file.path(config$data$figures, "01_pairwise_ANI_density_plot.pdf")
outfile_clusters <- file.path(config$basedir, "all_species_clusters.tsv")

# ---------------------- 1. Read and plot distances ----------------------
message("Reading mash distances from: ", mash_file)
mash_distances <- read.table(mash_file, sep="\t", header=FALSE)
colnames(mash_distances) <- c("Species", "Distance")

mash_distances <- mash_distances %>%
  dplyr::mutate(Similarity = 1 - Distance) %>%
  dplyr::mutate(Species = ifelse(grepl("Enterobacter", Species), "Enterobacter spp.", Species))

p <- ggplot(mash_distances, aes(x=Similarity, fill=Species, color=Species)) +
  geom_density(alpha = 0.1, adjust = 3) +
  coord_cartesian(xlim=c(0.95, 1)) +
  scale_fill_viridis_d(option = "C") +
  scale_color_viridis_d(option = "C") +
  labs(
    title = "Density of Pairwise MASH Distances by Species",
    x = "MASH Pairwise Distance",
    y = "Density"
  ) +
  theme_bw()

message("Saving plot to: ", density_plot_fp)
ggsave(density_plot_fp, p, width=6, height=4, units="in")

# ---------------------- 2. Cluster extraction ----------------------
message("Processing distance files in: ", distance_dir)
files <- list.files(distance_dir, pattern = "_dist_clean\\.txt$", full.names = TRUE)

all_clusters <- list()

for (f in files) {
  message("Processing: ", basename(f))
  species <- sub("_dist_clean.*", "", basename(f))
  df <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  if (nrow(df) == 0 || ncol(df) < 3) next
  colnames(df)[1:3] <- c("sample1", "sample2", "distance")
  
  samples <- unique(c(df$sample1, df$sample2))
  n <- length(samples)
  distmat <- matrix(NA, nrow=n, ncol=n, dimnames=list(samples, samples))
  diag(distmat) <- 0
  
  for (i in seq_len(nrow(df))) {
    s1 <- df$sample1[i]
    s2 <- df$sample2[i]
    d <- df$distance[i]
    distmat[s1, s2] <- d
    distmat[s2, s1] <- d
  }
  distmat[is.na(distmat)] <- max(df$distance, na.rm=TRUE) + 1
  
  distmat_dist <- as.dist(distmat)
  hc <- hclust(distmat_dist, method="complete")
  clust <- cutree(hc, h=threshold)
  outdf <- data.frame(
    sample = names(clust),
    cluster = as.integer(clust),
    species = species,
    stringsAsFactors = FALSE
  )
  
  keep <- names(which(table(clust) >= minsize))
  filtered <- outdf[outdf$cluster %in% keep, c("sample", "species", "cluster")]
  all_clusters[[species]] <- filtered
}

result <- do.call(rbind, all_clusters)
outfile_clusters <- file.path(distance_dir, "all_species_clusters.tsv")
message("Writing combined cluster file to: ", outfile_clusters)
write.table(result, file=outfile_clusters, sep="\t", row.names=FALSE, quote=FALSE)

message("✅ Done.")
