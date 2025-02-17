This directory compiles the annotation and tree files necessary to make the IS expansion phylogeny figures in Figure 2.

# several files were downloaded from iTol Templates

# loaded ../../IS_expansion_phylogeny/gtdbtk.bac120.decorated.tree

01 - Representative tree

All node labels
```bash
[is_expansion_tree] touch 01.enterococcus_tree_labels.txt
[is_expansion_tree] cat labels_template.txt itol_labels.txt > 01.enterococcus_tree_labels.txt # remove header line after DATA with vim
```

Pruned nodes
```bash
[is_expansion_tree] touch 01.enterococcus_tree_prune.txt
[is_expansion_tree] cat prune.txt representative_nodes.txt > 01.enterococcus_tree_prune.txt
```

Abundant IS family heatmap
```bash
[is_expansion_tree] touch 01.enterococcus_tree_heatmap.txt
[is_expansion_tree] cat dataset_heatmap_template.txt Abundant_IS_family_counts_per_reference.tsv > 01.enterococcus_tree_heatmap.txt # set MIN_COLOR and MAX_COLOR
# add COLOR_MAX #38761d COLOR_MIN #FFFFFF and FIELD_LABELS ISL3 IS30 IS256 IS3 IS6 IS110 to file with vim 
```

02 - E. faecium subclade tree

Pruned nodes
```bash
cat prune.txt efm_ISL3_prune_nodes.txt > 02.efm_clade_tree_prune.txt # remove header 'x' line from Data
```

Symbols
```bash
cat dataset_symbols_template.txt efm_clade_symbols.txt > efm_clade_tree_symbols.txt # convert tab to comma, remove header line from Data
```

ISL3 counts
```bash
cat dataset_simplebar_template.txt efm_ISL3_counts.txt > 02.efm_clade_tree_bar.txt # convert tab to comma, remove header line from Data
```
