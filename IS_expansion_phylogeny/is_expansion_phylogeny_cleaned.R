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
library(ggtreeExtra)
library(ggstar)
library(ggbeeswarm)
library(fuzzyjoin)

setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/IS_expansion_phylogeny")

# read representative genome csv, format for input into gtdb-tk ----------

rep_genomes <- read.delim("enterococci_genomes_sorted.csv", sep = ",")
input <- rep_genomes %>%
  dplyr::mutate(accession = gsub("\\./rep_genomes/", "", File.name)) %>%
  dplyr::mutate(superkingdom = "Bacteria") %>%
  dplyr::mutate(phylum = "Bacillota") %>%
  dplyr::mutate(class = "Bacilli") %>%
  dplyr::mutate(order = "Lactobacillales") %>%
  dplyr::mutate(family = "Enterococcaceae") %>%
  dplyr::mutate(ScientificName = paste(Genus, Species))

# Your existing data frame
df <- data.frame(
  ScientificName = c(
    "Lactococcus garvieae", "Carnobacterium maltaromaticum", 
    "Vagococcus lutrae", "Melissococcus plutonius", 
    "Enterococcus faecalis", "Enterococcus moraviensis", 
    "Enterococcus haemoperoxidus", "Enterococcus caccae", 
    "Enterococcus dispar", "Enterococcus asini", 
    "Enterococcus gallinarum", "Enterococcus casseliflavus", 
    "Enterococcus cecorum", "Enterococcus columbae", 
    "Enterococcus italicus", "Enterococcus sulfureus", 
    "Enterococcus saccharolyticus", "Tetragenococcus halophilus", 
    "Enterococcus pallens", "Enterococcus gilvus", 
    "Enterococcus avium", "Enterococcus raffinosus", 
    "Enterococcus malodoratus", "Enterococcus phoeniculicola", 
    "Enterococcus mundtii", "Enterococcus durans", 
    "Enterococcus hirae", "Enterococcus villorum", 
    "Enterococcus faecium"
  )
)

# Create a data frame for taxid mapping
taxid_mapping <- data.frame(
  ScientificName = c(
    "Lactococcus garvieae", "Carnobacterium maltaromaticum", 
    "Vagococcus lutrae", "Melissococcus plutonius", 
    "Enterococcus faecalis", "Enterococcus moraviensis", 
    "Enterococcus haemoperoxidus", "Enterococcus caccae", 
    "Enterococcus dispar", "Enterococcus asini", 
    "Enterococcus gallinarum", "Enterococcus casseliflavus", 
    "Enterococcus cecorum", "Enterococcus columbae", 
    "Enterococcus italicus", "Enterococcus sulfureus", 
    "Enterococcus saccharolyticus", "Tetragenococcus halophilus", 
    "Enterococcus pallens", "Enterococcus gilvus", 
    "Enterococcus avium", "Enterococcus raffinosus", 
    "Enterococcus malodoratus", "Enterococcus phoeniculicola", 
    "Enterococcus mundtii", "Enterococcus durans", 
    "Enterococcus hirae", "Enterococcus villorum", 
    "Enterococcus faecium"
  ),
  taxid = c(
    1363, 2751, 81947, 33970, 
    1351, 155617, 155618, 317735, 
    44009, 57732, 1353, 37734, 
    44008, 1355, 246144, 1356, 
    41997, 51669, 160454, 160453, 
    33945, 71452, 71451, 154621, 
    53346, 53345, 1354, 112904, 
    1352
  )
)

# Join the taxid to your existing data frame
df <- df %>%
  left_join(taxid_mapping, by = "ScientificName")

input <- input %>%
  dplyr::left_join(df, by = "ScientificName") %>%
  dplyr::select(accession, taxid, superkingdom, phylum, class, order, family, Genus, Species, Strain) %>%
  # rename Genus, Species, Strain to lowercase
  dplyr::rename(genus = Genus, species = Species, strain = Strain) %>%
  # remove .fna from accession's
  dplyr::mutate(accession = gsub("\\.fna", "", accession))

# format like this: genome_1 d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
gtdb_input <- input %>%
  # Create assembly_accession column from accession
  dplyr::mutate(assembly_accession = accession) %>%
  # Create taxonomy string with proper format and prefixes
  dplyr::mutate(taxonomy = paste0(
    "d__", superkingdom, ";p__", phylum, ";c__", class, ";o__", order, 
    ";f__", family, ";g__", genus, ";s__", species)) %>%
  # Select only the relevant columns
  dplyr::select(assembly_accession, taxonomy)

# You can then write gtdb_input to a file if needed:
write.table(gtdb_input, "gtdb_input_cleaned.tsv", row.names = FALSE, quote = FALSE, sep = "\t")


# Gilmore Enterococcus tree -----------------------------------------------

# Load the tree and input genomes
tree <- read.tree("enterococcus_genomes_gtdb_v2/gtdbtk.bac120.decorated.tree")
ref_genomes <- read.delim("gtdb_input_cleaned.tsv", sep = "\t")

# prepare iTol templates --------------------------------------------------

# extract leaf labels from tree
tip_nodes <- tree$tip.label

# prepare iTol labels
itol_labels <- data.frame(NODE_ID = tip_nodes)

# Improve formatting of labels for better visualization 
itol_labels <- itol_labels %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID %in% ref_genomes$assembly_accession, NODE_ID, NA)) %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID == "Enterococcus_faecalis_OG1RF_genomic", "Enterococcus faecalis OG1RF", LABEL)) %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID == "Enterococcus_faecalis_V583_genomic", "Enterococcus faecalis V583", LABEL)) %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID == "e_mundtii", "Enterococcus mundtii", LABEL)) %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID == "GCF_000250945.1_ASM25094v1_genomic", "Enterococcus faecium AUS0004", LABEL)) %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID == "AHXU01", "Enterococcus faecium EnGen0043", LABEL)) %>%
  dplyr::mutate(LABEL = ifelse(NODE_ID == "AHWP01", "Enterococcus faecium EnGen0015", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("GCF_[0-9]+\\.[0-9]+_", "", LABEL)) %>% # remove GCF_ and version number
  dplyr::mutate(LABEL = gsub("_genomic", "", LABEL)) %>% # remove _genomic
  dplyr::mutate(LABEL = gsub("_[^_]*$", "", LABEL)) %>% # strip the last underscore and everything after it
  dplyr::mutate(LABEL = gsub("Ente_", "Enterococcus ", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("Tetra_", "Tetragenococcus ", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("Vago_", "Vagococcus ", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("Carnob_", "Carnobacterium ", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("Lacto_", "Lactococcus ", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("Melissococcus_", "Melissococcus ", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("_[^_]*$", "", LABEL)) %>% # strip the last underscore and everything after it
  dplyr::mutate(LABEL = gsub("disp", "dispar", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("cass", "casseliflavus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("gall", "gallinarum", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("ceco", "cecorum", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("raffi", "raffinosus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("colu", "columbae", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("pall", "pallens", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("malod", "malodoratus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("ital", "italicus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("vill", "villorum", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("asin", "asini", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("phoe", "phoeniculicola", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("sacc", "saccharolyticus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("haem", "haemoperoxidus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub(" cocc", " caccae", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("gilv", "gilvus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("mora", "moraviensis", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("sulf", "sulfureus", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("garv", "garvieae", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("malt", "maltaromaticum", LABEL)) %>%
  dplyr::mutate(LABEL = gsub("halo", "halophilus", LABEL)) %>%
  dplyr::mutate(LABEL = ifelse(is.na(LABEL), NODE_ID, LABEL)) %>%
  dplyr::filter(LABEL != "EF013", LABEL != "EF_B_72") # these two have distances from E. faecalis > .15

# add CLASS information for the representative genomes
# using BAPS clusters from Lebreton et al 2017, Table S1
cluster_out <- c("Vagococcus lutrae", "Carnobacterium maltaromaticum", "Lactococcus garvieae")
cluster_1 <- c("Melissococcus plutonius")
cluster_2 <- c("Enterococcus faecalis OG1RF", "Enterococcus faecalis V583", "Enterococcus moraviensis", 
               "Enterococcus haemoperoxidus", "Enterococcus caccae")
cluster_3 <- c("Enterococcus dispar", "Enterococcus gallinarum", "Enterococcus casseliflavus", "Enterococcus cecorum", "Enterococcus columbae",
               "Enterococcus italicus", "Enterococcus saccharolyticus", 
               "Enterococcus sulfureus", "Tetragenococcus halophilus")
cluster_4 <- c("Enterococcus asini")
cluster_5 <- c("Enterococcus gilvus", "Enterococcus pallens", "Enterococcus raffinosus",
               "Enterococcus avium", "Enterococcus malodoratus")
cluster_6 <- c("Enterococcus phoeniculicola")
cluster_7 <- c("Enterococcus villorum", "Enterococcus hirae", "Enterococcus durans",  "Enterococcus mundtii", 
               "Enterococcus faecium AUS0004", "Enterococcus faecium EnGen0043", "Enterococcus faecium EnGen0015")

itol_labels <- itol_labels %>%
  dplyr::mutate(CLASS = case_when(
    LABEL %in% cluster_out ~ "outgroup",
    LABEL %in% cluster_1 ~ "one",
    LABEL %in% cluster_2 ~ "two",
    LABEL %in% cluster_3 ~ "three",
    LABEL %in% cluster_4 ~ "four",
    LABEL %in% cluster_5 ~ "five",
    LABEL %in% cluster_6 ~ "six",
    LABEL %in% cluster_7 ~ "seven",
    TRUE ~ NA
  ))

# write itol_labels to file
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/is_expansion_tree/")
write.table(itol_labels, file = "itol_labels.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# using tree distances, get clusters to each representative ---------------
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/IS_expansion_phylogeny")
dist_matrix <- read.delim("gtdbtk.bac120.decorated.tree_v2.dist", sep = ",", row.names = 1, check.names = F)
rownames(dist_matrix) <- gsub(" ", "_", rownames(dist_matrix))
colnames(dist_matrix) <- gsub(" ", "_", colnames(dist_matrix))

# get references
representative_genomes <- read.delim("isescan/representative_genomes.txt", sep = "\t", header = F, check.names = F)
colnames(representative_genomes) <- c("Sample") 
representative_genomes$Sample <- gsub(".fna", "", representative_genomes$Sample)

# subset distance matrix to references
dist_matrix <- dist_matrix[ , representative_genomes$Sample]

# Convert the distance matrix to long format
dist_long <- dist_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Genome") %>%
  pivot_longer(-Genome, names_to = "Reference", values_to = "Distance")

# filter out EF013 and EF_B_72
dist_long <- dist_long %>%
  dplyr::filter(Genome != "EF013", Genome != "EF_B_72")

# for each genome, get the reference it is closest to
closest_references <- dist_long %>%
  dplyr::group_by(Genome) %>%
  slice_min(Distance) %>%  # Get the row with the minimum distance
  ungroup()

# Extract the closest reference names
closest_references <- closest_references %>%
  dplyr::select(Genome, Reference) %>%  # Keep only necessary columns
  dplyr::mutate(Assembly = str_extract(Genome, "^[^.]+\\.[0-9]+")) %>%
  dplyr::mutate(Assembly = ifelse(is.na(Assembly), Genome, Assembly)) %>%
  dplyr::mutate(Assembly = gsub("\\.v1", "\\.1", Assembly)) %>%  # Standardize the assembly names
  dplyr::mutate(Assembly = gsub("\\.v2", "\\.2", Assembly)) %>%  # Replace underscores with spaces
  dplyr::left_join(itol_labels, by = c("Reference" = "NODE_ID"))  # Add the class information

# plot number of genomes closest to each reference
p <- ggplot(closest_references, aes(x = LABEL)) +
  geom_bar() +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of genomes closest to each reference",
       x = "Representative", y = "Genomes") +
  geom_text(stat = 'count', aes(label = ..count..),  # Add text labels for counts
            size = 4, vjust = -0.5)
p
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/is_expansion_tree")
ggsave("Number_of_genomes_by_closest_reference.svg", p, width = 10, height = 5, units = "in")

# read in ISEScan results, get IS counts by family for each genome --------
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/IS_expansion_phylogeny/isescan")
# Load and preprocess data
representative_genome_ises <- read.delim("all_isescan_sum.tsv", sep = "\t", header = TRUE)
representative_genome_ises$Sample <- gsub(".fna", "", representative_genome_ises$Sample)

representative_genome_ises_per_family <- representative_genome_ises %>%
  dplyr::filter(IS_Family != "total") %>%  # Remove the total row
  dplyr::left_join(itol_labels, by = c("Sample" = "NODE_ID")) %>% # Add class information
  dplyr::group_by(Sample, IS_Family, LABEL) %>%
  dplyr::summarize(count = sum(IS_count), .groups = "drop")

# Define the desired order of LABEL values to match tree
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
family_order <- c(
  "ISL3", "IS256", "IS30", "IS6", "IS3", "IS110", "IS982", "IS1182", "IS200/IS605", 
  "IS91", "IS1380", "IS701", "IS4", "IS66", "IS21", "IS630", "IS1634", "new",
  "IS5", "IS607"
)

# Ensure IS_Family and LABEL exist before calling complete()
representative_genome_ises_per_family <- representative_genome_ises_per_family %>%
  dplyr::mutate(
    IS_Family = factor(IS_Family, levels = family_order),
    LABEL = factor(LABEL, levels = label_order)
  ) %>%
  dplyr::rename(sample = Sample)

# Fill in missing LABEL-IS_Family combinations with count = 0
representative_genome_ises_per_family <- representative_genome_ises_per_family %>%
  tidyr::complete(LABEL, IS_Family, fill = list(count = 0))

# Create heatmap
p <- ggplot(representative_genome_ises_per_family, aes(y = LABEL, x = IS_Family, fill = count)) +
  geom_tile() +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IS Family Counts by Genome",
       x = "IS Family", y = "Genomes") +
  scale_fill_viridis_c()

# Print the plot
p

# print heatmap values to tsv, with IS families as columns and genomes as rows, need to pivot wider
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/figures/is_expansion_tree")
representative_genome_ises_per_family_wide <- pivot_wider(representative_genome_ises_per_family, 
                                                          names_from = IS_Family, 
                                                          values_from = count) %>%
  dplyr::select(-LABEL) %>%
  dplyr::select(sample, all_of(family_order))

# Save the wide format data
write.table(representative_genome_ises_per_family_wide, "IS_family_counts_by_genome.tsv", sep = " ", 
            quote = FALSE, row.names = FALSE)

# repeat this but now for all results
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/IS_expansion_phylogeny")
all_ises <- read.delim("isescan_per_sample.csv", sep = ",", header = T)

# filter to the Enterococcus set
# to do so, must convert .1, .2, .3, .4, etc. to _1 and so forth
all_ises_renamed <- all_ises %>%
  dplyr::mutate(sample = gsub("\\.(\\d+)", "_\\1", sample))

# bind representative genome IS data to all IS data
reps_to_bind <- representative_genome_ises_per_family %>%
  dplyr::select(sample, IS_Family, count) %>%
  dplyr::rename(family = IS_Family)

all_ises_renamed <- all_ises_renamed %>%
  dplyr::bind_rows(reps_to_bind)

all_ises_per_family <- closest_references %>%
  dplyr::left_join(all_ises_renamed, by = c("Genome" = "sample"))
  
# verify correct number of genomes
length(unique(all_ises_per_family$Genome)) #2134 

plot_data <- all_ises_per_family %>%
  dplyr::group_by(LABEL, family) %>%
  dplyr::summarize(average = as.numeric(mean(count)),
                   stdev = as.numeric(sd(count)),
                   max_count = as.numeric(max(count)),
                   max_genome = Genome[which.max(count)],
                   n = as.numeric(n())) %>%
  dplyr::left_join(itol_labels, by = "LABEL") %>%
  dplyr::mutate(LABEL = factor(LABEL, levels = label_order)) %>%
  dplyr::mutate(family = factor(family, levels = family_order)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(family != "NA")
  
# Create the heatmap
p <- ggplot(plot_data, aes(y = LABEL, x = family, fill = average)) +
  geom_tile() +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IS Family Counts by Genome",
       x = "IS Family", y = "Genomes") +
  scale_fill_viridis_c()

# Print the plot
p

# Create a new dataframe with all combinations of Reference and IS_Family
all_combinations <- expand.grid(LABEL = label_order, family = family_order)

# Merge the original data with the complete combinations, filling NAs with 0, relinking the NODE_ID
IS_counts_complete <- all_combinations %>%
  dplyr::left_join(plot_data, by = c("LABEL", "family")) %>%
  dplyr::select(-NODE_ID, -CLASS) %>%
  dplyr::left_join(itol_labels, by = "LABEL") %>%
  dplyr::mutate(
    average = ifelse(is.na(average), 0, average),
    stdev = ifelse(is.na(stdev), 0, stdev),
    n = ifelse(is.na(n), 0, n)
  )

# pivot IS_counts_complete wider
all_genome_ises_per_family_wide <- IS_counts_complete %>%
  dplyr::select(NODE_ID, family, average) %>%
  dplyr::mutate(average = as.numeric(average)) %>%
  dplyr::filter(!is.na(NODE_ID)) %>%
  pivot_wider(names_from = family, 
              values_from = average)

# Save the wide format data
write.table(all_genome_ises_per_family_wide, "IS_family_counts_per_reference_all.tsv", sep = " ", 
            quote = FALSE, row.names = FALSE)

# repeat but for only abundant IS families
# list of abundant IS families for plotting
abundant_is_families <- c("ISL3", "IS30", "IS256", "IS3", "IS6", "IS110")

all_genome_ises_abundant_family_wide <- IS_counts_complete %>%
  dplyr::select(NODE_ID, family, average) %>%
  dplyr::mutate(average = as.numeric(average)) %>%
  dplyr::filter(!is.na(NODE_ID)) %>%
  pivot_wider(names_from = family, 
              values_from = average) %>%
  dplyr::mutate(other = rowSums(select(., -NODE_ID, -all_of(abundant_is_families[-length(abundant_is_families)])), na.rm = TRUE)) %>%
  dplyr::select(NODE_ID, any_of(abundant_is_families))

# Save the wide format data
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/")
write.table(all_genome_ises_abundant_family_wide, "figures/is_expansion_tree/Abundant_IS_family_counts_per_reference.tsv", sep = " ", 
            quote = FALSE, row.names = FALSE)

# write closest references and ISE counts per genome and reference to files ----------

# write closest references to file
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/IS_expansion_phylogeny/")
write.table(closest_references, "all_enterococcus_closest_references.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# write ISE counts per reference to file
write.table(IS_counts_per_family_per_reference, "all_enterococcus_IS_counts_per_reference.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# write all E. faecalis genome names to a file -----------------------------------
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/IS_expansion_phylogeny/")
efs_genomes <- closest_references %>%
  dplyr::filter(grepl("faecalis", LABEL)) %>%
  dplyr::select(Genome) %>%
  dplyr::mutate(Genome = paste0(Genome, ".fna"))

write.table(efs_genomes, "efs_genomes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# write all E. faecium genome names to a file -----------------------------

efm_genomes <- closest_references %>%
  dplyr::filter(grepl("faecium", LABEL)) %>%
  dplyr::select(Genome) %>%
  dplyr::mutate(Genome = paste0(Genome, ".fna"))

write.table(efm_genomes, "efm_genomes.txt", sep = "\t", quote = FALSE, row.names = FALSE)



# get clade ranges for clade B, clade A1, clade A2 ------------------------

efm_ranges <- closest_references %>%
  dplyr::filter(grepl("faecium", LABEL)) %>%
  dplyr::mutate(Clade = case_when(
    grepl("aus0004", LABEL, ignore.case = T) ~ "clade A1",
    grepl("engen0043", LABEL, ignore.case = T) ~ "clade A2",
    grepl("engen0015", LABEL, ignore.case = T) ~ "clade B",
    TRUE ~ "other"
  )) %>%
  dplyr::select(Genome, Clade)

# color palette
clade_colors <- moma.colors(n=3, "Picasso")

# symbols based on clade
efm_clade_symbols <- efm_ranges %>%
  dplyr::mutate(feature = 2) %>% # 2 is circle
  dplyr::mutate(symbol_color = case_when(Clade == "clade A1" ~ clade_colors[1],
                                         Clade == "clade A2" ~ clade_colors[2],
                                         Clade == "clade B" ~ clade_colors[3],
                                         TRUE ~ "#000000")) %>%
  dplyr::mutate(symbol_size = 1) %>%
  dplyr::mutate(symbol_fill = 1) %>%
  dplyr::mutate(symbol_position = 1) %>%
  dplyr::mutate(label = "efm clade") %>%
  dplyr::select(Genome, feature, symbol_size, symbol_color, symbol_fill, symbol_position, label)

# write to file
setwd("/Volumes/lab_asbhatt/mpgriesh/projects/transmission/data/sra/efm/")
write.table(efm_clade_symbols, "figures/is_expansion_tree/efm_clade_symbols.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE)


# write ISL3 counts for Efm genomes to file -------------------------------

isl3_counts <- efm_ranges %>%
  dplyr::left_join(
    all_ises_per_family %>%
      dplyr::filter(family == "ISL3") %>%
      dplyr::select(Genome, count),  # Keep only relevant columns
    by = "Genome"
  ) %>%
  dplyr::mutate(count = tidyr::replace_na(count, 0)) %>%  # Replace missing counts with 0
  dplyr::select(Genome, count)  # Keep only required columns

write.table(isl3_counts, "figures/is_expansion_tree/efm_ISL3_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# prune nodes for those with ISL3 counts
isl3_prune_nodes <- isl3_counts %>%
  dplyr::select(Genome) %>%
  write.table("figures/is_expansion_tree/efm_ISL3_prune_nodes.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# clade specific mean/sd/genomes ------------------------------------------

a1 <- efm_ranges %>%
  dplyr::left_join(all_ises_per_family %>% dplyr::filter(family == "ISL3"), by = "Genome") %>%
  dplyr::mutate(count = tidyr::replace_na(count, 0)) %>%
  dplyr::group_by(Clade) %>%
  dplyr::summarize(mean = mean(count), sd = sd(count), n = n())
