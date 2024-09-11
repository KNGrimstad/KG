<h1>
 KG - A Collection of scRNA-seq Functions<a href="https://github.com/KNGrimstad/KG"><img align="right" src="images/KG_logotype.png" width="148.59" height="165.4575"/></a>
</h1>


<!-- badges: start -->
[![Static Badge](https://img.shields.io/badge/Version-0.1.0-lightblue)](https://github.com/KNGrimstad/KG/releases/tag/v0.1.0)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/KNGrimstad/KG)](https://github.com/KNGrimstad/KG)
[![GitHub License](https://img.shields.io/github/license/KNGrimstad/KG)](https://github.com/KNGrimstad/KG?tab=MIT-1-ov-file)
[![Static Badge](https://img.shields.io/badge/repo%20status-active-lightgreen)](https://www.repostatus.org/#active)
<!-- badges: end -->

This repository is composed of functions and wrappers designed to make life easier for the owner of the repo, and maybe for you, too. The functions included are all, in one way or another, related to processing and analysis of scRNA-seq data. 

Some of these functions are exceptionally base, while others are used to construct figures of varying degrees of complexity and customizability. Note that these functions were created with and for use with the Seurat package prior to the v5 update. Functions may be updated in the future, but Seurat v5 is currently not supported.

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

These functions are either original plotting functions or wrappers for existing ones, aiming to generate more clean and publication-friendly figures. <img align="right" src="images/example_plots4.png" width = "265" height = "340">

| **Function** | **Description** |
| --- | --- |
| **`KG_dimplot`** | Make a clean dimplot for puclications. |
| **`KG_heatmap`** | Fancy heatmap for cell-by-cell or cluster-average expression of genes. |
| **`KG_isotype_pie`** | Pie chart of Ig isotype distributions in clusters. |
| **`KG_subtype_pie`** | Same as above, but for Ig subtypes |
| **`KG_3DPlot`** | A classic dimension plot in 3D. |
| **`KG_dotplot`** | A basic dot plot, albeit subjectively nicer. |
| **`KG_clone_tree`** | Wrapper to make clone trees from V(D)J-seq data. |
| **`KG_percent_cells`** | Plots percentage of cells in each cluster. |
| **`KG_isotype_dimplot`** | Dimension plots that highlight cells by their IgH isotype. |

## Data wrangling

The functions in this section aim to make life a bit easier for working with scRNA-seq data, in particular focusing on B-cell datasets. 

| **Function** | **Description** |
| --- | --- |
| **`KG_remove_BCR`** | Removes the BCR genes from the dataset, see [**Sundell et al.(2022)**](https://doi.org/10.1093/bfgp/elac044). |
| **`KG_import_seurat_clusters`** | Imports Seurat annotations to a monocle3 dataset. |
| **`KG_filter_pbmcs`** | Plots expression for common PBMC genes and prompts you which clusters you want to keep. |
| **`KG_gene_to_isotype`** | Converts C calls from V(D)J-seq data to Ig isotypes. |
| **`KG_isotype_switch`** | Converts Ig isotype to switch status (unswitched/switched). |
| **`KG_get_trajectory`** | Extracts the coordinates from monocle3-generated trajectories. |
| **`KG_filter_vdj`** | Filter V(D)J-seq data frames for downstream BCR analyses. |
| **`KG_matchGEXtoVDJ`** | Matched V(D)J-seq data to paired GEX data. |

## Helpful wrappers

The functions here are all wrappers to existing functions, including ones in Seurat and the Immcantation framework. 

| **Function** | **Description** |
| --- | --- |
| **`KG_clones`** | Wrapper for calling clones in V(D)J-seq data. Wraps Immcantation framework functions. |
| **`KG_cluster_multi´** | Wrapper that runs Seurat::FindClusters for a range of resolutions, and prompts you to select the one you want to go ahead with. |
| **`KG_integrate_CCA`** | Wrapper for Seurat's CCA integration workflow. | 
| **`KH_hc`** | Wrapper for Seurat's hierarchical clustering workflow. |
| **`KG_project`** | Wrapper for Seurat's functions for projecting datasets onto a reference in dimensionally reduced space. | 
| **`KG_project_3D`** | Same as above, but generates a 3D plot. |
| **`KG_prep_clones`** | Streamlines base clonotype processing from the Immcantation framework. |
| **`KG_multi_gene_plot`** | Plots gene expression (feature plot or violin plot) for a list of genes, one by one. |

## Misc. helpers

The functions in this category are very basic and only strive to reduce the time spent writing the same code over and over again. Perfect for lazy days. 

| **Function** | **Description** |
| --- | --- |
| **`KG_add_suffix`** | Adds a suffix to all strings supplied. Nice for cell barcodes. |
| **`KG_remove_suffix`** | Same as above, but for removing a suffix. |
| **`KG_fetch_genes`** | Returns all genes in a Seurat object with the supplied gene name prefix. Nice for pesky HLA genes. | 
| **`KG_plot_to_pdf`** | Writes plots to PDF, with customizable layouts. |
| **`KG_percent_isotypes`** | Calculates the percent of Ig isotypes in each Seurat identity. |

# Repository Activity
![Alt](https://repobeats.axiom.co/api/embed/9baa88e7488279b7170d442f240ed5cc46abfd5a.svg "Repobeats analytics image")
