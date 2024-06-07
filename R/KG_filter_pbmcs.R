# Filter data based on cell type markers.

#' A function for differentiating B cells, T cells, and monocytes in PBMC scRNA-seq data.
#' @param seurat_object A Seurat object.
#' @param features The genes to use for cell type determination.
#' @param scale.factor Scale factor to use for normalization following filtering.
#' @export
#' @examples
#' KG_filter_pbmcs(B_cell_dataset)
KG_filter_pbmcs = function(seurat_object,
                           features = c("CD19", "CD3E", "CD14"),
                           scale.factor = 10000){

  suppressPackageStartupMessages({
    require(Seurat)
    require(patchwork)
    require(base)})

  plots = list()
  plots[["Dimplot"]] = Seurat::DimPlot(seurat_object, pt.size = 2, label = T, label.size = 6)
  plots[["VlnPlots"]] = suppressWarnings(Seurat::VlnPlot(seurat_object, features = features,
                                                         ncol = 1))

  print(patchwork::wrap_plots(plots, ncol = 2))

  keep_clusters = readline(prompt = "Which clusters do you want to keep? (separated only by comma): ")
  cat("Subsetting and reprocessing ")
  keep_clusters = strsplit(keep_clusters, ',')
  keep_clusters = as.numeric(unlist(keep_clusters))
  seurat_object = base::subset(seurat_object, idents = keep_clusters)
  Seurat::DefaultAssay(seurat_object) = "RNA"
  seurat_object = Seurat::NormalizeData(object = seurat_object,
                                        scale.factor = scale.factor,
                                        verbose = F)
  seurat_object = Seurat::ScaleData(object = seurat_object,
                                    features = rownames(seurat_object),
                                    verbose = F)
  return(seurat_object)
}
