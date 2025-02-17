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
library(ggnetwork)
library(ggraph)


setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/lr_meta_sv/')

# plotting theme ----------------------------------------------------------

theme_custom <- theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme_bw() +
  theme(
    text = element_text(color = "black", face = "plain"),
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(color = "black", size = 7),
    plot.title = element_text(color = "black", size = 7),
    legend.position = "none"
  )
# long read metaflye qc and extract efm contigs ---------------------------
setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/lr_meta/metaflye')
lr_info <- read.delim('all_assembly_info.txt', sep = '\t')
circ_contigs <- lr_info %>%
  dplyr::filter(circular == 'Y')

long_circ_contigs <- circ_contigs %>%
  dplyr::filter(length > 1000000)

# write Sample,seq_name to tab delimited file
long_circ_contigs %>%
  dplyr::select(Sample, seq_name) %>%
  write.table("long_circ_contigs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# read in ska distances ---------------------------------------------------
setwd('/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/lr_meta_sv/')
distances <- read.delim('05.SKA/chromosome_results/distances.tsv', sep = '\t')
clusters <- read.delim('05.SKA/chromosome_results/clusters.tsv', sep = '\t')
df <- distances %>%
  dplyr::mutate(sample_1 = str_split_i(Sample.1, "_contig", 1),
                sample_2 = str_split_i(Sample.2, "_contig", 1)) %>%
  dplyr::mutate(sample_1 = gsub("_chromosome", "", sample_1),
                sample_2 = gsub("_chromosome", "", sample_2)) %>%
  dplyr::mutate(contig_1 = str_split_i(Sample.1, "_contig_", 2),
                contig_2 = str_split_i(Sample.2, "_contig_", 2)) %>%
  dplyr::mutate(contig_1 = gsub("_multi.fasta", "", contig_1),
                contig_2 = gsub("_multi.fasta", "", contig_2)) %>%
  dplyr::select(-Sample.1, -Sample.2) %>%
  dplyr::mutate(comp_1 = paste0(sample_1, "_", contig_1),
                comp_2 = paste0(sample_2, "_", contig_2)) %>%
  dplyr::mutate(Patient_1 = str_split_i(sample_1, "_", 1),
                Patient_2 = str_split_i(sample_2, "_", 1)) %>%
  dplyr::filter(Matches > 0)
  
p <- df %>%
  dplyr::filter(Patient_1 == Patient_2) %>%
  dplyr::filter(Matches >= 2000000) %>%
  ggplot(aes(x = comp_1, y = comp_2, fill = Jaccard.Index)) +
  geom_tile() +
  facet_grid(Patient_1 ~ Patient_2, scales = "free") +
  scale_fill_viridis_c() +
  theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) +
  theme_custom
p
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics")
ggsave("lr_meta_ska_per_patient_jaccard.pdf", plot = p, width = 8, height = 6, units="in")

# get absolute abundance data for each sample ------------------------------
setwd('/Users/matthewgrieshop/OneDrive - Stanford/MSTP/Bhatt/Transmission/LR_metagenomics')
sampled <- read.delim("efm_filtered.txt", sep = "\t")

abs_abundance <- read.delim("efm_abs_abundance.csv", sep = ",")
abs_abundance <- abs_abundance %>%
  dplyr::mutate(p5id = str_split_i(SampleID, "_", 1),
                date = as.Date(str_split_i(fullID, "_", 2)),
                total_efm = e_faecium_rel_abun * total_16S_rRNA_copies_per_wet_g_stool) %>%
  dplyr::mutate(Sequenced = ifelse(SampleID %in% sampled$SampleID, T, F)) %>%
  dplyr::group_by(p5id) %>%
  dplyr::mutate(
    first_sample_date = min(date, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(days_since_first = as.numeric(date - first_sample_date))

# filter for p5id with at least 2 samples where e_faecium_rel_abun >= 0.1
filtered <- abs_abundance %>%
  dplyr::filter(p5id != "") %>%
  dplyr::filter(e_faecium_rel_abun >= 0.1) %>%
  dplyr::group_by(p5id) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(p5id)
  
# plot relative abundance over time per patient
p <- abs_abundance %>%
  dplyr::filter(p5id %in% filtered$p5id) %>%
  ggplot(aes(x = days_since_first, y = e_faecium_rel_abun)) +
  geom_line() +
  geom_point(aes(x = days_since_first, y = e_faecium_rel_abun, color = Sequenced, size=5)) +
  facet_wrap(~p5id, scales = "free_x") +
  scale_color_manual(values = moma.colors("Smith", 2)) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  labs(x = "Days since first sample", y = "Relative abundance of E. faecium") +
  theme_custom
p

setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics")
# save pdf
ggsave("sampled_efm_rel_abun_over_time.pdf", plot = p, width = 8, height = 8, units="in")


# plot relative abundance over time colored by patient --------------------

p <- abs_abundance %>%
  dplyr::filter(p5id %in% filtered$p5id) %>%
  ggplot(aes(x = days_since_first, y = e_faecium_rel_abun)) +
  geom_line(aes(color = p5id)) +
  geom_point(aes(x = days_since_first, y = e_faecium_rel_abun)) +
  scale_color_manual(values = moma.colors("Smith", 18)) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  labs(x = "Days since first sample", y = "Relative abundance of E. faecium") +
  theme_custom
p

# save pdf
ggsave("all_sampled_efm_rel_abun_over_time.pdf", plot = p, width = 8, height = 8, units="in")

# collect assembly stats and sourmash stats to generate plot ---------------------------------
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/lr_meta_sv")

nanostats <- read.delim("07.Stats/all_nanostats.tsv", sep = "\t")
nanostats <- nanostats %>%
  dplyr::mutate(Sample = gsub("nanoplot_", "", Sample)) %>%
  dplyr::mutate(Sample = gsub("_post", "", Sample)) %>%
  dplyr::mutate(patient = str_split_i(Sample, "_", 1))

median((nanostats$n50/1000)) # 9.181 K
median(nanostats$Median_Qual) # 23.1

# plot read_N50 vs median read quality
p <- nanostats %>%
  ggplot(aes(x = n50/1000, y = Median_Qual, color = patient)) +
  geom_point() +
  scale_color_manual(values = moma.colors("Smith", 16)) +
  theme_bw() +
  geom_vline(data = nanostats, aes(xintercept = median(n50)/1000), linetype = "dashed", color = "red") +
  ylim(0, 25) +
  xlim(0, 15) +
  labs(x = "Read N50", y = "Median read quality") +
  theme_custom
p

# save pdf
ggsave("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics/nanoplot_n50_vs_median_qual.pdf", plot = p, width = 4, height = 4, units="cm")

stats <- read.delim("07.Stats/all_assembly_info.txt", sep = "\t")

sourmash <- read.delim("01.Sourmash/results.singletons.txt", sep = ",")
colnames(sourmash) <- c("sample","status","proportion","cANI","lineage","name")

# sourmash, split sample into patient and contig
sourmash <- sourmash %>%
  dplyr::mutate(p5id = str_split_i(sample, "_contig", 1),
                contig = paste0("contig", str_split_i(sample, "_contig", 2))) %>%
  dplyr::mutate(Patient = str_split_i(p5id, "_", 1)) %>%
  dplyr::left_join(stats, by = c("p5id" = "Sample", "contig" = "seq_name"))

enterococcus <- sourmash %>%
  dplyr::filter(grepl("Enterococcus", lineage, ignore.case = T)) 

# group by name, plot Coverage vs Length
p <- enterococcus %>%
  ggplot(aes(y = coverage, x = length/1000000, color = name)) +
  geom_point() +
  scale_color_manual(values = moma.colors("Smith", 9)) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) +
  labs(y = "Coverage", x = "Contig Length (Mb)") +
  theme_custom
p

# bar plot of p5id with >2.5Mb,faecium contig vs. not
summary_df <- enterococcus %>%
  dplyr::group_by(p5id) %>%
  dplyr::summarise(has_faecium = any(name == "faecium" & length > 2500000)) %>%
  dplyr::ungroup() %>%
  dplyr::count(has_faecium)

# Create the bar plot
p <- ggplot(summary_df, aes(x = has_faecium, y = n, fill = has_faecium)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = c("FALSE" = "No", "TRUE" = "Yes")) +
  labs(x = "Has faecium contig > 2.5Mb", y = "Number of p5id values", title = "Presence of faecium contigs > 2.5Mb") +
  theme_bw() +
  scale_fill_manual(values = moma.colors(n=2,"Liu")) +
  theme_custom

# save pdf
ggsave("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics/closed_efaecium_contig_summary.pdf", plot = p, width = 4, height = 6, units="cm")

# plot coverage vs length for all contigs, coloring faecium contigs
p <- sourmash %>%
  dplyr::mutate(plot_color = ifelse(grepl("faecium", name, ignore.case = T), "faecium", "other")) %>%
  dplyr::mutate(coverage = ifelse (coverage > 2500, 2500, coverage)) %>% # cap coverage at 2500
  ggplot(aes(y = coverage, x = length/1000000, color = plot_color)) +
  geom_point() +
  scale_color_manual(values = moma.colors("Smith", 2)) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) +
  labs(y = "Coverage", x = "Contig Length (Mb)") +
  theme_custom
p

# save pdf
ggsave("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics/all_contigs_coverage_vs_length.pdf", plot = p, width = 8, height = 5, units="cm")

# which samples don't?
faecium_status <- enterococcus %>%
  dplyr::group_by(p5id) %>%
  dplyr::mutate(has_faecium = any(name == "faecium" & length > 2500000)) %>%
  dplyr::select(p5id, has_faecium) %>%
  dplyr::distinct() %>%
  dplyr::mutate(patient = str_split_i(p5id, "_", 1)) %>%
  dplyr::group_by(patient) %>%
  dplyr::summarize(count = sum(has_faecium)) %>%
  dplyr::filter(count > 1) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::mutate(patient = factor(patient))

# join with treatment information
faecium_status_meta <- faecium_status %>%
  dplyr::left_join(abs_abundance, by = c("patient" = "p5id"))

# plot SNP distance from SKA within each patient --------------------------

efm_chromosomes <- enterococcus %>%
  dplyr::filter(length > 2500000 & name == "faecium") %>%
  dplyr::mutate(contig = gsub("contig_", "", contig)) %>%
  dplyr::mutate(key = paste(p5id, contig, sep = "_"))

efm_chromosomes %>%
  dplyr::select(p5id, contig) %>%
  dplyr::mutate(key = paste(p5id, "contig", contig, sep = "_")) %>%
  dplyr::select(key) %>%
  write.table("08.FaeciumChromosomes/efm_chromosomes.txt", quote = F, row.names = F, col.names = F)

# filter df to only include comparisons of efm_chromosomes (p5id and contig match in both columns)
efm_chromosome_df <- df %>%
  dplyr::filter(comp_1 %in% efm_chromosomes$key & comp_2 %in% efm_chromosomes$key)

# filter to cases where Patient_1 == Patient_2
efm_chromosome_df_within_patient <- efm_chromosome_df %>%
  dplyr::filter(Patient_1 == Patient_2)

# get first sample per patient
first_sample <- abs_abundance %>%
  dplyr::filter(Sequenced) %>%
  dplyr::filter(p5id %in% efm_chromosomes$Patient) %>%
  dplyr::group_by(p5id) %>%
  dplyr::summarize(first_sample = min(SampleID))

# join first_sample to efm_chromosome_df_within_patient
efm_chromosome_df_within_patient <- efm_chromosome_df_within_patient %>%
  dplyr::left_join(first_sample, by = c("Patient_1" = "p5id"))

# plot dot plot of SNP distance index per patient for E. faecium contigs
p <- efm_chromosome_df_within_patient %>%
  dplyr::filter(sample_1 == first_sample) %>%
  dplyr::mutate(similarity = 100*(1-SNP.distance)) %>%
  ggplot(aes(x = first_sample, y = similarity, color = Patient_2)) +
  geom_point(size = 5) +
  scale_color_manual(values = moma.colors("Smith", 12)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "None") +
  labs(y = "% Identity", x = "First sample sequenced") +
  coord_flip() +
  ylim(99.995, 100) +
  geom_hline(yintercept = 99.999, linetype = "dashed", color = "red") +
  theme_custom
p

# save as pdf
ggsave("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics/efm_chromosome_similarity_within_patient.pdf", plot = p, width = 5, height = 5, units="cm")


# visualize .dot graph with ggnetwork -------------------------------------

# plot cluster graph
chromosome_edge <- read.csv("05.SKA/chromosome_results/chromosome_results.edgelist.csv")
chromosome_network <- ggnetwork(chromosome_edge)

# filter out any 'repeat' values
chromosome_network <- chromosome_network %>%
  dplyr::mutate(Patient = str_split_i(vertex.names, "_chromosome", 1)) %>%
  dplyr::mutate(Contig = str_split_i(vertex.names, "_chromosome_", 2)) %>%
  dplyr::mutate(vertex.names = gsub("_chromosome", "", vertex.names)) %>%
  dplyr::mutate(key = paste(Patient, Contig, sep = "_")) %>%
  dplyr::mutate(key = gsub("_contig", "", key))

# clean vertex.names
chromosome_network$vertex.names <- gsub("_contig_[0-9]+$", "", chromosome_network$vertex.names)

# get only efm chromosomes
efm_chromosome_keys <- unique(c(efm_chromosome_df$comp_1, efm_chromosome_df$comp_2))

# join metadata and filter to only include efm chromosomes
chromosome_network_meta <- chromosome_network %>%
  dplyr::left_join(abs_abundance, by = c("vertex.names" = "SampleID")) %>%
  dplyr::filter(key %in% efm_chromosome_keys)

# get connected nodes
chromosome_edge_count <- chromosome_network %>%
  dplyr::group_by(x, y, vertex.names) %>%
  dplyr::summarize(edge_count = n())

cluster_labels <- chromosome_network_meta %>%
  distinct(x, y, vertex.names, .keep_all = TRUE)

# Create the plot
ggplot(chromosome_network_meta) +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend, alpha = weight, linewidth = weight), color = "gray90") +  # Draw edges
  geom_nodes(aes(x = x, y = y, fill = label), size = 5) +  # Draw nodes
  geom_text_repel(data = cluster_labels, aes(x = x, y = y, label = vertex.names)) +  # Apply filtering only for labels
  scale_fill_manual(values = moma.colors('Smith', n = 18)) +
  theme_bw()

# noticed that some of the edges go to nowhere, this is because I filtered for Efm chromosomes. To fix, rerunning SKA on just the Efm chromosomes and plotting that

# visualize .dot of Efm chromosomes alone with ggnetwork ------------------

# first get the distances/clusters from ska output
distances <- read.delim('05.SKA/faecium_results/results.distances.tsv', sep = '\t')
clusters <- read.delim('05.SKA/faecium_results/results.clusters.tsv', sep = '\t')
df <- distances %>%
  dplyr::mutate(sample_1 = str_split_i(Sample.1, "_contig", 1),
                sample_2 = str_split_i(Sample.2, "_contig", 1)) %>%
  dplyr::mutate(sample_1 = gsub("_chromosome", "", sample_1),
                sample_2 = gsub("_chromosome", "", sample_2)) %>%
  dplyr::mutate(contig_1 = str_split_i(Sample.1, "_contig_", 2),
                contig_2 = str_split_i(Sample.2, "_contig_", 2)) %>%
  dplyr::mutate(contig_1 = gsub("_multi.fasta", "", contig_1),
                contig_2 = gsub("_multi.fasta", "", contig_2)) %>%
  dplyr::select(-Sample.1, -Sample.2) %>%
  dplyr::mutate(comp_1 = paste0(sample_1, "_", contig_1),
                comp_2 = paste0(sample_2, "_", contig_2)) %>%
  dplyr::mutate(Patient_1 = str_split_i(sample_1, "_", 1),
                Patient_2 = str_split_i(sample_2, "_", 1)) %>%
  dplyr::filter(Matches > 0)

# plot cluster graph
chromosome_edge <- read.csv("05.SKA/faecium_results/results.edgelist.csv")
chromosome_network <- ggnetwork(chromosome_edge)

# filter out any 'repeat' values
chromosome_network <- chromosome_network %>%
  dplyr::mutate(Patient = str_split_i(vertex.names, "_chromosome", 1)) %>%
  dplyr::mutate(Contig = str_split_i(vertex.names, "_chromosome_", 2)) %>%
  dplyr::mutate(vertex.names = gsub("_chromosome", "", vertex.names))

# clean vertex.names
chromosome_network$vertex.names <- gsub("_contig_[0-9]+$", "", chromosome_network$vertex.names)

# get connected nodes
chromosome_edge_count <- chromosome_network %>%
  dplyr::group_by(x, y, vertex.names) %>%
  dplyr::summarize(edge_count = n())

cluster_labels <- chromosome_network_meta %>%
  distinct(x, y, vertex.names, .keep_all = TRUE)

# Create the plot
ggplot(chromosome_network_meta) +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend, alpha = weight, linewidth = weight), color = "gray90") +  # Draw edges
  geom_nodes(aes(x = x, y = y, fill = label), size = 5) +  # Draw nodes
  geom_text_repel(data = cluster_labels, aes(x = x, y = y, label = vertex.names)) +  # Apply filtering only for labels
  scale_fill_manual(values = moma.colors('Smith', n = 18)) +
  theme_bw()

# visualize 1-SNP.distance as heatmap -------------------------------------

# first join cluster information so they can be used for faceting
clusters <- clusters %>%
  dplyr::mutate(sample = str_split_i(ID, "_contig", 1)) %>%
  dplyr::mutate(sample = gsub("_chromosome", "", sample)) %>%
  dplyr::mutate(contig = str_split_i(ID, "_contig_", 2)) %>%
  dplyr::mutate(contig = gsub("_multi.fasta", "", contig))

df <- df %>%
  dplyr::left_join(clusters, by = c("sample_1" = "sample", "contig_1" = "contig")) %>%
  dplyr::left_join(clusters, by = c("sample_2" = "sample", "contig_2" = "contig"), suffix = c(".1", ".2"))

p <- df %>%
  dplyr::filter(Cluster__autocolour.1 == Cluster__autocolour.2) %>%
  ggplot(aes(x = comp_1, y = comp_2, fill = SNPs)) +
  geom_tile() +
  facet_wrap(~Cluster__autocolour.1, scales = "free") +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) +
  theme_custom +
  theme(legend.position = "bottom")
p

# save as pdf
ggsave("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/lr_metagenomics/all_ska_clusters.pdf", plot = p, width = 15, height = 15, units="cm")


