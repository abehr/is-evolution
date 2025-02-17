# Packages and directory ----------------------------------------------------------------

library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(reshape2)
library(forcats)
library(viridis)
library(ggthemes)
library(ggmsa)
library(msa)
library(ggtree)
library(treeio)
library(ggrepel)
library(phangorn)
library(ape)
library(Biostrings)
library(tinytex)
library(MoMAColors)
library(ggridges)
library(pheatmap)
library(readr)
library(ggnewscale)
library(tidytree)
library(TreeTools)
library(ggtreeExtra)
library(ggstar)
library(ggbeeswarm)
library(fuzzyjoin)

setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/poppunk")

# efaecalis poppunk analysis --------------------------------------------------------
setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/pangraph/efaecalis/')

# get poppunk cluster information per Sample
clusters <- read.delim("clusters/clusters_clusters.csv", sep = ",", header = T)

# get MLST output to attempt to assign ST per cluster, looking for biggest cluster
mlst <- read.delim("mlst.txt", sep = "\t", header = F)
colnames(mlst) <- c("Sample", "Scheme", "ST", "gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL")
mlst <- mlst %>%
  ## correct Sample name to reflect editing done by poppunk (converts '.' to '_' )
  dplyr::mutate(Sample = gsub(".fna", "", Sample)) %>%
  dplyr::mutate(Sample = gsub("inputs/", "", Sample)) %>%
  dplyr::mutate(Sample = gsub("\\.", "_", Sample))

# join clusters and mlst
clusters_mlst <- mlst %>%
  dplyr::left_join(clusters, by = c("Sample" = "Taxon"))

# get cluster sizes, with most common ST assigned
cluster_sizes <- clusters_mlst %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n(), ST = names(which.max(table(ST))))

# plot cluster sizes (10 or more) and label with ST
p <- cluster_sizes %>%
  dplyr::filter(n > 9) %>%
  dplyr::arrange(desc(n)) %>%
  ggplot(aes(x = reorder(Cluster, n), y = n)) +
  geom_bar(stat = "identity") +
  theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Cluster Sizes",
       x = "Cluster", y = "Number of Genomes") +
  geom_text(aes(label = paste0("ST", ST), vjust = -0.5, size = 5))
p

# quick look at the tree
nwk <- read.tree("viz/viz_core_NJ.nwk")
p <- ggtree(nwk, layout = 'circular') + theme_tree2() + geom_text2(aes(label=node), size = 2)
p

df <- as_tibble(nwk)

# prepare iTol templates --------------------------------------------------

# extract leaf labels from tree
tip_nodes <- nwk$tip.label

# prepare iTol labels
itol_labels <- data.frame(NODE_ID = tip_nodes) 

# add CLASS information equal to Cluster
itol_labels <- itol_labels %>%
  dplyr::left_join(clusters_mlst, by = c("NODE_ID" = "Sample")) %>%
  dplyr::mutate(LABEL = gsub("Enterococcus_faecalis_", "", NODE_ID)) %>%
  dplyr::mutate(LABEL = ifelse(str_detect(LABEL, "_ASM"),
                               str_extract(LABEL, ".*(?=_ASM)"),
                               LABEL)) %>%
  dplyr::mutate(LABEL = str_replace(LABEL, "(?<!_B)_1", "\\.1")) %>%
  dplyr::mutate(LABEL = str_replace(LABEL, "(?<!_B)_2", "\\.2")) %>%
  dplyr::mutate(LABEL = str_replace(LABEL, "(?<!_B)_3", "\\.3")) %>%
  dplyr::select(NODE_ID, LABEL, Cluster) %>%
  dplyr::rename(CLASS = Cluster)

# change to figures/poppunk/efaecalis
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/poppunk/efaecalis")
# write itol_labels to file
write.table(itol_labels, file = "itol_labels.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# find nodes defining the range for cluster 1_113, 2, 8

colored_clades <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% c("1_113", "2", "8")) %>%
  dplyr::left_join(cluster_sizes, by = c("CLASS" = "Cluster")) %>%
  dplyr::select(CLASS, ST) %>%
  dplyr::distinct() %>%
  dplyr::mutate(min_node = case_when(CLASS == "2" ~ "EF_B_65",
                                     CLASS == "8" ~ "EF_B_15",
                                     CLASS == "1_113" ~ "EF_B_49",
                                     TRUE ~ "NA")) %>%
  dplyr::mutate(max_node = case_when(CLASS == "2" ~ "GCF_036279155_1",
                                     CLASS == "8" ~ "GCF_022212865_1",
                                     CLASS == "1_113" ~ "GCF_029278545_1",
                                     TRUE ~ "NA")) %>%
  dplyr::mutate(type = "clade") %>%
  dplyr::mutate(style = "normal") %>%
  dplyr::mutate(scale = 2) %>%
  dplyr::mutate(ST = paste0("ST", ST)) %>%
  dplyr::mutate(color_1 = moma.colors(n=3, "Liu")) %>%
  dplyr::mutate(color_2 = moma.colors(n=3, "Liu")) %>%
  dplyr::mutate(line_color = "#000000") %>%
  dplyr::mutate(line_style = "solid") %>%
  dplyr::mutate(line_scale = "1") %>%
  dplyr::mutate(label_color = "#000000") %>%
  dplyr::mutate(label_scale = "1") %>%
  dplyr::mutate(label_style = "normal") %>%
  dplyr::select(min_node, max_node, color_1, color_2, line_color, line_style, line_scale, ST, label_color, label_scale, label_style)

# write to file
write.table(colored_clades, file = "colored_clades.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# symbols for SHC samples
shc_symbols <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(grepl("FM|EF", NODE_ID)) %>%
  dplyr::left_join(cluster_sizes, by = c("CLASS" = "Cluster")) %>%
  dplyr::select(NODE_ID, CLASS, ST) %>%
  dplyr::mutate(feature = 2) %>% # 2 is circle
  dplyr::mutate(symbol_color = "#ff0000") %>%
  dplyr::mutate(symbol_size = 3) %>%
  dplyr::mutate(symbol_fill = 1) %>%
  dplyr::mutate(symbol_position = 1) %>%
  dplyr::mutate(label = "SHC") %>%
  dplyr::select(NODE_ID, feature, symbol_size, symbol_color, symbol_fill, symbol_position, label)

# write to file
write.table(shc_symbols, file = "shc_symbols.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# print lists of NODE_ID to files, save in working directory --------------
setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/pangraph/efaecalis/')

cluster_1 <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% c("8")) %>%
  dplyr::select(LABEL) %>%
  dplyr::mutate(LABEL = paste0(LABEL, ".fna"))

# write to file
write.table(cluster_1, file = "cluster_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cluster_2 <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% c("2")) %>%
  dplyr::select(LABEL) %>%
  dplyr::mutate(LABEL = paste0(LABEL, ".fna"))

# write to file
write.table(cluster_2, file = "cluster_2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cluster_3 <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% c("1_113")) %>%
  dplyr::select(LABEL) %>%
  dplyr::mutate(LABEL = ifelse(grepl("V583", LABEL), "Enterococcus_faecalis_V583_genomic", LABEL)) %>%
  dplyr::mutate(LABEL = paste0(LABEL, ".fna"))

# write to file
write.table(cluster_3, file = "cluster_3.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)





# efaecium poppunk analysis -----------------------------------------------
setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/pangraph/efaecium/')

# get poppunk cluster information per Sample
clusters <- read.delim("clusters/clusters_clusters.csv", sep = ",", header = T)

# get MLST output to attempt to assign ST per cluster, looking for biggest cluster
mlst <- read.delim("mlst.txt", sep = "\t", header = F)
colnames(mlst) <- c("Sample", "Scheme", "ST", "gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL")
mlst <- mlst %>%
  ## correct Sample name to reflect editing done by poppunk (converts '.' to '_' )
  dplyr::mutate(Sample = gsub(".fna", "", Sample)) %>%
  dplyr::mutate(Sample = gsub("inputs/", "", Sample)) %>%
  dplyr::mutate(Sample = gsub("\\.", "_", Sample))

# join clusters and mlst
clusters_mlst <- mlst %>%
  dplyr::left_join(clusters, by = c("Sample" = "Taxon"))

# get cluster sizes, with most common ST assigned
cluster_sizes <- clusters_mlst %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n(), ST = names(which.max(table(ST)))) %>%
  dplyr::filter(n > 9) %>%
  dplyr::filter(!is.na(Cluster))

# plot cluster sizes (10 or more) and label with ST
p <- cluster_sizes %>%
  dplyr::arrange(desc(n)) %>%
  ggplot(aes(x = reorder(Cluster, n), y = n)) +
  geom_bar(stat = "identity") +
  theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Cluster Sizes",
       x = "Cluster", y = "Number of Genomes") +
  geom_text(aes(label = paste0("ST", ST), vjust = -0.5, size = 5))
p

# quick look at the tree
nwk <- read.tree("viz/viz_core_NJ.nwk")
p <- ggtree(nwk, layout = 'circular') + theme_tree2() + geom_text2(aes(label=node), size = 2)
p

df <- as_tibble(nwk)

df <- df %>%
  dplyr::left_join(clusters_mlst, by = c("label" = "Sample")) %>%
  dplyr::mutate(color = ifelse(Cluster %in% c("1_21", "39", "40", "41", "42", "43"), Cluster, NA),
                size = ifelse(Cluster %in% c("1_21", "39", "40", "41", "42", "43"), 1, 0.2))

tree_data <- as.phylo(df)
p <- ggtree(tree_data, layout = 'circular', aes(color = color, size = size)) %<+% df +
  geom_rootedge(0.001) +
  #geom_text(aes(label = label), hjust = 1, vjust = 1, size = 2) +
  scale_color_manual(values = MoMAColors::moma.colors("vonHeyl", n = 7)) +
  guides(color = guide_legend(title = NULL)) +  # Remove legend title +
  theme_tree2(legend.position = "top")
p

# prepare iTol templates --------------------------------------------------

# extract leaf labels from tree
tip_nodes <- nwk$tip.label

# prepare iTol labels
itol_labels <- data.frame(NODE_ID = tip_nodes) 

# add CLASS information equal to Cluster
itol_labels <- itol_labels %>%
  dplyr::left_join(clusters_mlst, by = c("NODE_ID" = "Sample")) %>%
  dplyr::mutate(LABEL = gsub("Enterococcus_faecium_", "", NODE_ID)) %>%
  dplyr::mutate(LABEL = ifelse(str_detect(LABEL, "_ASM"),
                               str_extract(LABEL, ".*(?=_ASM)"),
                               LABEL)) %>%
  dplyr::mutate(LABEL = str_replace(LABEL, "(?<!_B)_1", "\\.1")) %>%
  dplyr::mutate(LABEL = str_replace(LABEL, "(?<!_B)_2", "\\.2")) %>%
  dplyr::mutate(LABEL = str_replace(LABEL, "(?<!_B)_3", "\\.3")) %>%
  dplyr::select(NODE_ID, LABEL, Cluster) %>%
  dplyr::rename(CLASS = Cluster)

# change to figures/poppunk/efaecium
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/poppunk/efaecium")
# write itol_labels to file
write.table(itol_labels, file = "itol_labels.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

colored_clades <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% c("40")) %>%
  dplyr::select(CLASS, ST) %>%
  dplyr::distinct() %>%
  dplyr::mutate(min_node = "FM053") %>%
  dplyr::mutate(max_node = "FM_B_15") %>%
  dplyr::mutate(type = "clade") %>%
  dplyr::mutate(style = "normal") %>%
  dplyr::mutate(scale = 2) %>%
  dplyr::mutate(color_1 = moma.colors(n=1, "Liu")) %>%
  dplyr::mutate(color_2 = moma.colors(n=1, "Liu")) %>%
  dplyr::mutate(line_color = "#000000") %>%
  dplyr::mutate(line_style = "solid") %>%
  dplyr::mutate(line_scale = "1") %>%
  dplyr::mutate(label_color = "#000000") %>%
  dplyr::mutate(label_scale = "1") %>%
  dplyr::mutate(label_style = "normal") %>%
  dplyr::select(min_node, max_node, color_1, color_2, line_color, line_style, line_scale, CLASS, label_color, label_scale, label_style)

# write to file
write.table(colored_clades, file = "colored_clades.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# color samples in clusters
# Define color mapping
class_colors <- setNames(moma.colors(n = 2, "Liu"), c("39", "40"))

colored_labels <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% names(class_colors)) %>%
  dplyr::mutate(feature = "label",
                type = "node",
                label_color = class_colors[CLASS],
                background_color = class_colors[CLASS],
                label_scale = "1",
                label_style = "normal") %>%
  dplyr::select(NODE_ID, feature, type, label_color, label_scale, label_style, background_color)

# write to file
write.table(colored_labels, file = "colored_labels.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# color SHC samples
shc_labels <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(grepl("FM|EF", NODE_ID)) %>%
  dplyr::select(NODE_ID, CLASS, ST) %>%
  dplyr::mutate(feature = "label") %>%
  dplyr::mutate(type = "node") %>%
  dplyr::mutate(label_color = "#000000") %>%
  dplyr::mutate(label_scale = "1") %>%
  dplyr::mutate(label_style = "normal") %>%
  dplyr::mutate(background_color = "#ff0000") %>%
  dplyr::select(NODE_ID, feature, type, label_color, label_scale, label_style, background_color)

# write to file
write.table(shc_labels, file = "colored_shc_labels.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# symbols for SHC samples
shc_symbols <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(grepl("FM|EF", NODE_ID)) %>%
  dplyr::select(NODE_ID, CLASS, ST) %>%
  dplyr::mutate(feature = 2) %>% # 2 is circle
  dplyr::mutate(symbol_color = "#ff0000") %>%
  dplyr::mutate(symbol_size = 3) %>%
  dplyr::mutate(symbol_fill = 1) %>%
  dplyr::mutate(symbol_position = 1) %>%
  dplyr::mutate(label = "SHC") %>%
  dplyr::select(NODE_ID, feature, symbol_size, symbol_color, symbol_fill, symbol_position, label)

# write to file
write.table(shc_symbols, file = "shc_symbols.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# print lists of NODE_ID to files, save in working directory --------------
setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/pangraph/efaecium/')

cluster_1 <- itol_labels %>%
  dplyr::left_join(df, by = c("NODE_ID" = "label")) %>%
  dplyr::filter(CLASS %in% c("40")) %>%
  dplyr::filter(grepl("FM|EF", NODE_ID)) %>%
  dplyr::select(LABEL) %>%
  dplyr::mutate(LABEL = paste0(LABEL, ".fna"))

# write to file
write.table(cluster_1, file = "cluster_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
