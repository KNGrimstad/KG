<h1>
  KG - A Collection of Random scRNA-Seq Functions&nbsp;<img align = "right" src = "images/KG_logotype.png" width = "114.3" height = "127.275">
</h1>

This repository is composed of functions designed to make life easier for the owner of the repo. The functions included are all, in one way or another, related to processing and analysis of scRNA-seq data. 

Further descriptions will come.

---
## Installation
Install the most recent version of this package: 
```r
devtools::install("KNGrimstad/KG")
```
---
# Functions Currently Included in this Package
## Fancy plotting functions
KG_heatmap - make a fancy heatmap for cell-by-cell or cluster-average expression of genes.<img align = "right" src = "images/example_plots.png" width = "114" height = "115">

KG_isotype_pie - make

KG_subtype_pie

KG_3DPlot

KG_clone_tree

## Data wrangling
KG_remove_BCR

KG_import_seurat_clusters

KG_filter_bcr

KG_gene_to_isotype

KH_isotype_to_switch

## Helpful(?) wrappers
KG_clones

KG_cluster_multi

KG_integrate_CCA



## Minor helpers
KG_add_suffix

KG_remove_suffix

KG_fetch_genes

KG_plot_to_pdf
