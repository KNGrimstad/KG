<h1>
  KG - A Collection of Random scRNA-Seq Functions&nbsp;<img align = "right" src = "images/KG_logotype.png" width = "114.3" height = "127.275">
</h1>

<!-- badges: start -->
![GitHub R package version](https://img.shields.io/github/r-package/v/KNGrimstad/KG)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/KNGrimstad/KG)
![GitHub License](https://img.shields.io/github/license/KNGrimstad/KG)
![Static Badge](https://img.shields.io/badge/repo%20status-active-lightgreen)
<!-- badges: end -->
This repository is composed of functions and wrappers designed to make life easier for the owner of the repo, and maybe for you, too. The functions included are all, in one way or another, related to processing and analysis of scRNA-seq data. 

Some of these functions are exceptionally base, while others are used to construct figures of varying degrees of complexity and customizability. 

If you are working with lymphocyte scRNA-seq data, consider giving our [*protocol paper*](https://doi.org/10.1093/bfgp/elac044) a read.

---
## Installation
Install the most recent version of this package: 
```r
devtools::install_github("KNGrimstad/KG")
```
---
# Functions Currently Included in this Package
## Fancy plotting functions 
<img align = "right" src = "images/example_plots4.png" width = "265" height = "340">

* **KG_dimplot** - Make clean dimplots for publications.

* **KG_heatmap** - Fancy heatmap for cell-by-cell or cluster-average expression of genes.

* **KG_isotype_pie** - Pie chart of Ig isotype distributions in clusters.

* **KG_subtype_pie** - Same as above, but for Ig subtypes.

* **KG_3DPlot** - A classic dimension plot in 3D

* **KG_dotplot** - Seurat's dot plot, but subjectively nicer.

* **KG_clone_tree** - Wrapper to make clone trees from V(D)J-seq data.

* **KG_percent_cells** - Plots percentage of cells in each cluster.

* **KG_isotype_dimplot** - Dimplots that highlight cells by their IGH isotype.

## Data wrangling
* **KG_remove_BCR** - Removes BCR genes from the dataset. 

* **KG_import_seurat_clusters** - Imports Seurat annotations to monocle3 object.

* **KG_filter_pbmcs** - Plots expression for common PBMC genes so you can filter out the cells you want.

* **KG_gene_to_isotype** - Converts c calls to main Ig isotypes.

* **KG_isotype_switch** - Converts Ig isotype to switch status (unswitched/switched).

* **KG_get_trajectory** - Exracts the coordinates from monocle3-generated trajectories.

* **KG_filter_vdj** - Filter AIRR format files from V(D)J-seq for downstream BCR analyses.

* **KG_matchGEXtoVDJ** - Matches V(D)J-seq data to paired GEX data.

## Helpful wrappers
* **KG_clones** - Wrapper to call clones in V(D)J-seq data.

* **KG_cluster_multi** - Wrapper to run Seurat's FindClusters for a range of resolutions.

* **KG_integrate_CCA** - Wrapper for Seurat's CCA integration workflow. 

* **KG_hc** - Wrapper for Seurat's hierarchical clustering.

* **KG_project** - Wrapper for Seurat's functions for projecting datasets onto a reference.

* **KG_project_3D** - Same as above, but generates a 3D dimension plot. 

* **KG_prep_clones** - Streamlines base clonotype processing from the Immcantation framework.

* **KG_multi_gene_plot** - Plots gene expression (feature or violoin) for a list of genes one by one.

## Misc. helpers
* **KG_add_suffix** - Simply adds a suffix to all strings supplied (optimal for cell barcodes).

* **KG_remove_suffix** - Same as above, but for removing a suffix.

* **KG_fetch_genes** - Returns all genes in a Seurat object with the supplied gene name prefix.

* **KG_plot_to_pdf** - Writes plots to PDF, with customizable layouts.

* **KG_percent_isotypes** - Calculates the percent of Ig isotypes in each Seurat identity.

