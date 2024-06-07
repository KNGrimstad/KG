#' Evaluate cluster resolutions with clustree
#'
#' This function runs the function FindClusters from the Seurat package, for a range of resolution values, plots the clusters for each resolution, and prompts the user to select an appropriate resolution, and returns the clustered Seurat object for the selected resolution.
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

  suppressPackageStartupMessages({
    require(ggtree)
    require(clustree)
    require(Seurat)
    require(SeuratObject)
    require(progress)})

  if(verbose){
    pb = progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta",
                                    total = length(resolutions),
                                    show_after = 0,
                                    complete = "+",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
  }

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
  if(verbose){
    cat("Finding neighbors\n")
  }
  seurat_object = Seurat::FindNeighbors(seurat_object,
                                        assay = assay,
                                        reduction = reduction,
                                        dims = dims,
                                        verbose = FALSE)

  resolutions = seq.int(min.res, max.res, increment)
  if(verbose){
    cat("Finding clusters\n")
  }
  for (i in resolutions){
    if(verbose){
      pb$tick()
    }
    seurat_object = Seurat::FindClusters(seurat_object, resolution = i, verbose = FALSE)
  }
  if(verbose){
    cat("Constructing cluster tree\n")
  }
  t = clustree::clustree(seurat_object@meta.data, prefix = paste(assay, "_snn_res.", sep = ""))
  plot(t)
  select_resolution = readline(prompt = "Select cluster resolution: ")

  if(verbose){
    cat(paste("Finding clusters at a resolution of ", select_resolution, "\n", sep = ""))
  }
  seurat_object = Seurat::FindClusters(seurat_object, resolution = as.numeric(select_resolution), verbose = FALSE)

  if(verbose){
    cat("Done\n")
  }
  return(seurat_object)
}
