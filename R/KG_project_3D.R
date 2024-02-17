#' Project and plot in 3D
#'
#' This function projects a Seurat object onto the clusters and UMAP coordinates from a reference, and plots the projected data in 3D.
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
#' @param pt.size Controls the size of plotted points/cells.
#' @param cols Determines the color palette to use for cells. If NULL, uses default colors.
#' @param bg_col Controls the background color of the plot.
#' @param min.dist Controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimize more accurately the local structure.
#' @param spread Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.
#' @export
#' @examples
#' # KG_project_3D(B_cell_dataset, reference_object)
KG_project_3D = function(seurat_object,
                         reference,
                         assay = "RNA",
                         reference.assay = "RNA",
                         reduction = "pca",
                         dims = 1:35,
                         min.dist = 0.3,
                         spread = 1,
                         reference.reduction = "pca",
                         reduction.model = "umap",
                         reference.idents = NULL,
                         pt.size = 5,
                         cols = NULL,
                         bg_col = "white"){
  require(Seurat)
  require(plotly)

  temp = seurat_object
  Seurat::DefaultAssay(temp) = assay

  # Prep reference data for 3D
  ref = Seurat::RunUMAP(reference,
                        dims = dims,
                        return.model = T,
                        n.components = 3,
                        min.dist = min.dist,
                        spread = spread)

  # Project data
  anchors = Seurat::FindTransferAnchors(reference = ref,
                                        query = temp,
                                        dims = dims,
                                        reference.reduction = reference.reduction)

  predictions = Seurat::TransferData(anchorset = anchors,
                                     refdata = ref[[reference_idents, drop = TRUE]],
                                     dims = dims)
  temp = Seurat::AddMetaData(temp, metadata = predictions)
  temp = Seurat::MapQuery(anchorset = anchors, reference = ref,
                          query = temp,
                          refdata = list(celltype = reference_idents),
                          reference.reduction = reduction,
                          reduction.model = "umap")
  Seurat::Idents(temp) = temp$predicted.celltype

  # Plot 3D
  df = data.frame(umap1 = temp@reductions$ref.umap@cell.embeddings[,
                                                                   1], umap2 = temp@reductions$ref.umap@cell.embeddings[, 2],
                  umap3 = temp@reductions$ref.umap@cell.embeddings[, 3], cell = Idents(temp))

  # Select color palette
  if(is.null(cols)){
    cols = hue_pal()(length(levels(df$cell)))
  }

  # Plot in 3D
  plotly::plot_ly(df, x = ~umap1,
                  y = ~umap2,
                  z = ~umap3,
                  color = ~cell,
                  marker = list(size = pt.size),
                  colors = cols) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = "UMAP_1"),
                        yaxis = list(title = "UMAP_2"),
                        zaxis = list(title = "UMAP_3"))) %>%
    layout(plot_bgcolor = bg_col) %>%
    layout(paper_bgcolor = bg_col)
}
