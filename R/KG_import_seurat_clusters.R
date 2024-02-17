#' Import Seurat-defined clusters to Monocle object
#' This function imports the clustering and dimension plot coordinates from a Seurat object and adds them to a Monocle3 object.
#' @param monocle_object A monocle object.
#' @param seurat_object A Seurat object.
#' @param dims The number of principal components to use.
#' @param plot_3d Logical, if TRUE, computes and plots in three dimensions.
#' @param assay Name of the assay from the Seurat object to use, "RNA" by default.
#' @param sparse Logical, if TRUE, attempts to increase distance between cells in the dimension plot.
#' @param min.dist If sparse = TRUE, controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimizes more accurately the local structure.
#' @param spread Effective scale of embedded points. If sparse = TRUE, determines, along with min.dist, how clumped embedded points are.
#' @param ident_col The identity column to use from the Seurat object. By default, 'seurat_clusters' is used.
#' KG_import_seurat_clusters(monocle_object, B_cell_dataset)
KG_import_seurat_clusters = function(monocle_object,
                                     seurat_object,
                                     dims = 35,
                                     plot_3d = FALSE,
                                     assay = "RNA",
                                     sparse = FALSE,
                                     min.dist = 0.3,
                                     spread = 1,
                                     ident_col = "seurat_clusters"){
  #DefaultAssay(seurat_object) = assay

  # For plotting in 3D
  if(plot_3d == TRUE){

    # Prep the Seurat object for 3D plotting
    require(plotly)
    require(scales)
    if(sparse == TRUE){
      message("Calculating 3D components")
      temp = RunUMAP(object = seurat_object, dims = 1:dims,
                     n.components = 3,
                     spread = spread,
                     min.dist = min.dist,
                     verbose = F,
                     assay = assay)
    } else {
      message("Calculating 2D components")
      temp = RunUMAP(object = seurat_object,
                     dims = 1:dims,
                     n.components = 3,
                     verbose = F,
                     assay = assay)
    }

  }

  # Create a base partitioning
  recreate.partition = c(rep(1, length(monocle_object@colData@rownames)))
  names(recreate.partition) = monocle_object@colData@rownames
  recreate.partition = as.factor(recreate.partition)
  monocle_object@clusters$UMAP$partitions = recreate.partition

  # Assign cluster information
  list_cluster = seurat_object[[ident_col, drop = TRUE]]
  monocle_object@clusters$UMAP$clusters = list_cluster

  # Assigning UMAP coordinate-cell embeddings
  if(plot_3d == TRUE){
    monocle_object@int_colData@listData$reducedDims$UMAP =
      temp@reductions$umap@cell.embeddings
  } else {
    monocle_object@int_colData@listData$reducedDims$UMAP =
      seurat_object@reductions$umap@cell.embeddings
  }


  return(monocle_object)
}
