#' Evaluate cluster resolutions with clustree
#'
#' This function runs the function FindClusters from the Seurat package, for a range of resolution values and returns a tree plot of the resulting clusters for each given resolution.
#' @param seurat_object A Seurat object.
#' @param max.res The maximum resolution for clustering.
#' @param min.res The minimum resolution for clustering.
#' @param increment A value determining by how much each increment in resolution should be.
#' @param plot If TRUE, returns a tree plot of resulting clusters for each resolution.
#' @param assay Name of the Seurat assay to use, "RNA" by default.
#' @param dims The number of dimensions to use for FindNeighbors.
#' @param reduction Name of the dimensionality reduction method to use, "pca" by default.
#' @export
#' @examples
#' KG_cluster_multi(B_cell_dataset)
KG_cluster_multi = function(seurat_object,
                            max.res = 1.0,
                            min.res = 0.1,
                            increment = 0.1,
                            plot = TRUE,
                            assay = "RNA",
                            dims = 1:35,
                            reduction = "pca"){
  require(ggtree)
  require(clustree)
  require(Seurat)
  require(SeuratObject)

  Seurat::DefaultAssay(seurat_object) = assay

  # Remove previous clustering if data is integrated
  if(nrow(unique(seurat_object[["orig.ident"]])) >1){
    cluster_columns = colnames(seurat_object@meta.data)[grepl(paste("^", assay, "_snn_res.", sep = ""),
                                                              colnames(seurat_object@meta.data))]

    for (col in cluster_columns){
      seurat_object@meta.data[[col]] = NULL
    }
    seurat_object = SeuratObject::UpdateSeuratObject(seurat_object)
  }
  # Perform clustering
  message("Finding neighbors")
  seurat_object = Seurat::FindNeighbors(seurat_object,
                                        assay = assay,
                                        reduction = reduction,
                                        dims = dims)

  resolutions = seq.int(min.res, max.res, increment)
  seurat_object = Seurat::FindNeighbors(seurat_object,
                                        dims = dims,
                                        assay = assay,
                                        reduction = reduction)
  for (i in resolutions){
    message(paste("Finding clusters for resolution", i, sep = " "))
    seurat_object = Seurat::FindClusters(seurat_object, resolution = i)
  }
  if(plot == TRUE){
      message("Constructing cluster tree")
      t = clustree::clustree(seurat_object@meta.data, prefix = paste(assay, "_snn_res.", sep = ""))
      plot(t)
  }
  return(seurat_object)
}
