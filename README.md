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
**KG_heatmap** - Fancy heatmap for cell-by-cell or cluster-average expression of genes.<img align = "right" src = "images/example_plots.png" width = "460" height = "460">

**KG_isotype_pie** - Pie chart of Ig isotype distributions in clusters.

**KG_subtype_pie** - Same as above, but for Ig subtypes.

**KG_3DPlot** - A classic dimension plot in 3.D

**KG_dotplot** - Seurat's dot plot, but subjectively nicer.

**KG_clone_tree** - Wrapper to make clone trees from V(D)J-seq data.

## Data wrangling
**KG_remove_BCR** - Removes BCR genes from the dataset. 

**KG_import_seurat_clusters** - Imports Seurat annotations to monocle3 object.

**KG_filter_bcr** - Filters an AIRR data frame from V(D)J-seq for quality.

**KG_gene_to_isotype** - Converts c calls to main Ig isotypes.

**KG_isotype_to_switch** - Converts Ig isotype to switch status (unswitched/switched).

## Helpful(?) wrappers
**KG_clones** - Wrapper to call clones in V(D)J-seq data.

**KG_cluster_multi** - Wrapper to run Seurat's FindClusters for a range of resolutions.

**KG_integrate_CCA** - Wrapper for Seurat's CCA integration workflow. 

## Minor helpers
**KG_add_suffix** - Simply adds a suffix to all strings supplied (optimal for cell barcodes).

**KG_remove_suffix** - Same as above, but for removing a suffix.

**KG_fetch_genes** - Returns all genes in a Seurat object with the supplied gene name prefix.

**KG_plot_to_pdf** - Writes plots to PDF, with customizable layouts.
