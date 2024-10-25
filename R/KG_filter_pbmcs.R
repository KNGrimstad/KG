# Filter data based on cell type markers.

#' A function for differentiating B cells, T cells, and monocytes in PBMC scRNA-seq data.
#' @param seurat_object A Seurat object.
#' @param features The genes to use for cell type determination.
#' @param dims Which dimensions to plot in the dim plot and potential feature plots.
#' @param ncol The number of columns for violin plots and feature plots.
#' @param plot Which type of plot to show the features as. Supported options are "violin" and "feature".
#' @param normalize Logical; whether data should be normalized after subsetting.
#' @param scale.factor Scale factor to use for normalization following filtering.
#' @param scale Logical; whether data should be scaled after subsetting.
#' @param find.var.feats Logical; whetehr new variable features should be identified after subsetting data.
#' @param n.var.feats How many variable features to identify.
#' @param runPCA Logical; whether Seurat::RunPCA should be applied after subsetting.
#' @param assay.reproc Name of the assay to use when re-processing.
#' @export
#' @examples
#' KG_filter_pbmcs(B_cell_dataset)
KG_filter_pbmcs_test = function(seurat_object,
                           features = c("CD19", "CD3E", "CD14"),
                           dims = c(1, 2),
                           assay = "RNA",
                           ncol = 1,
                           plot = "violin",
                           normalize = TRUE,
                           scale.factor = 10000,
                           scale = TRUE,
                           find.var.feats = FALSE,
                           n.var.feats = 2000,
                           runPCA = FALSE,
                           assay.reproc = "RNA"
                           ){

  suppressPackageStartupMessages({
    require(Seurat)
    require(patchwork)
    require(base)})

  plots = list()
  plots[["Dimplot"]] = KG_dimplot(seurat_object,
                                  pt.size = 2,
                                  stroke = 0.075,
                                  label = T,
                                  label.size = 7,
                                  dims = dims) + Seurat::NoLegend()
  if(plot == "violin"){
    plots[["VlnPlots"]] = suppressWarnings(Seurat::VlnPlot(seurat_object,
                                                           features = features,
                                                           ncol = ncol))
  } else if(plot == "feature"){
    plots[["Feature"]] = Seurat::FeaturePlot(seurat_object,
                                             features = features,
                                             ncol = ncol,
                                             cols = c("lightgrey", "red"),
                                             order = T,
                                             dims = dims)
  }

  print(patchwork::wrap_plots(plots, ncol = 2))

  # Select clusters to keep
  keep_clusters = readline(prompt = "Which clusters do you want to keep? (separated only by comma): ")
  cat("Subsetting and reprocessing \n")
  keep_clusters = strsplit(keep_clusters, ',')
  keep_clusters = as.numeric(unlist(keep_clusters))
  seurat_object = base::subset(seurat_object, idents = keep_clusters)
  Seurat::DefaultAssay(seurat_object) = assay.reproc

  # Normalize data?
  cat("Normalizing data\n")
  if(normalize){seurat_object = Seurat::NormalizeData(object = seurat_object,
                                                      scale.factor = scale.factor,
                                                      verbose = F)
  }
  # Scale data?
  if(scale){
    cat("Scaling data\n")
    seurat_object = Seurat::ScaleData(object = seurat_object,
                                      features = rownames(seurat_object),
                                      verbose = F)
  }
  # Find new variable features?
  if(find.var.feats){
    cat("Finding variable features\n")
    seurat_object = Seurat::FindVariableFeatures(seurat_object,
                                                 nfeatures = n.var.feats, verbose = FALSE)
  }
  # Run PCA?
  cat("Running PCA\n")
  if(runPCA){
    seurat_object = Seurat::RunPCA(seurat_object, verbose = FALSE)
  }
  cat("Done\n")
  return(seurat_object)
}
