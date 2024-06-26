#' Project data
#'
#' This function projects a Seurat object onto the clusters and UMAP coordinates from a reference.
#' @param seurat_object A Seurat object.
#' @param reference A Seurat object from which cluster identities and UMAP coordinates/embeddings should be pulled.
#' @param reduction Name of the dimensional reduction method to use as input, "pca" by default.
#' @param assay Name of the Seurat object assay to use.
#' @param reference.assay Name of the Seurat reference object assay to use.
#' @param dims The range of dimensions to use from linear dimensional reduction method.
#' @param reduction.model Name of dimensional reduction model to use for plotting, "umap" by default.
#' @param reference.reduction Name of dimensional reduction method to use as input from reference, "pca" by default.
#' @param reference.idents The identity variable to use from the reference for the projection. By default, uses the default identity of the reference object.
#' @param min.dist Controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimizes more accurately the local structure.
#' @param spread Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.
#' @export
#' @examples
#' KG_project(B_cell_dataset, B_cell_reference_dataset)
KG_project = function(seurat_object,
                      reference,
                      reduction = "pca",
                      assay = "RNA",
                      reference.assay = "RNA",
                      dims = 1:35,
                      min.dist = 0.3,
                      spread = 1,
                      reduction.model = "umap",
                      reference.reduction = "pca",
                      reference.idents = NULL){

  suppressPackageStartupMessages({
    require(Seurat)
    require(ggplot2)
  })

  temp = seurat_object
  Seurat::DefaultAssay(temp) = assay

  # Get reference model
  cat("Calculating components\n")
  ref = Seurat::RunUMAP(reference,
                        dims = dims,
                        return.model = T,
                        n.components = 2,
                        min.dist = min.dist,
                        spread = spread,
                        assay = reference.assay)
  # Project data
  cat("Finding transfer anchors\n")
  anchors = Seurat::FindTransferAnchors(reference = ref,
                                        query = temp,
                                        dims = dims,
                                        reference.reduction = reference.reduction,
                                        normalization.method = ifelse(assay == "SCT", "SCT", "LogNormalize"))

  cat("Transferring data\n")
  predictions = Seurat::TransferData(anchorset = anchors,
                                     refdata = ref[[reference.idents, drop = TRUE]],
                                     dims = dims)

  temp = Seurat::AddMetaData(temp, metadata = predictions)

  cat("Mapping query\n")
  temp = Seurat::MapQuery(anchorset = anchors,
                          reference = ref,
                          query = temp,
                          refdata = list(celltype = reference.idents),
                          reference.reduction = reduction,
                           reduction.model = "umap")
  Seurat::Idents(temp) = temp$predicted.celltype
  cat("Done\n")
  return(temp)
}
