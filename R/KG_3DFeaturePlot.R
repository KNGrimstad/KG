#' Feature plots, but in 3D!
#'
#' This function generates a feature plot in 3D.
#' @param seurat_object A Seurat object.
#' @param gene The gene for which the expression should be plotted.
#' @param reduction Name of the dimensional reduction method to use as input, "pca" by default.
#' @param dims The range of dimensions to use from linear dimensional reduction method.
#' @param assay Name of the Seurat object assay to use.
#' @param pt.size Controls the size of plotted points/cells.
#' @param stroke Size of border to data points.
#' @param cols Determines the colors to use for expression levels. Should be a vector of two colors.
#' @param bg.col Controls the background color of the plot.
#' @param min.dist Controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimize more accurately the local structure.
#' @param spread Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.
#' @param group.by Controls how cells should be grouped/colored. By default, uses the default identity of cells in the Seurat object.
#' @param trajectory_coords A vector containing coordinates for drawing trajectory lines.
#' @param trajectory_col What color to use for the trajectory.
#' @param lwd Line width for trajectory.
#' @export
#' @examples
#' KG_3DFeaturePlot(B_cell_dataset, gene = "CD27")
KG_3DFeaturePlot = function (seurat_object,
                             gene,
                             dims = 1:35,
                             reduction = "pca",
                             assay = "RNA",
                             pt.size = 7,
                             stroke = 1,
                             cols = c("lightgrey", "red"),
                             bg_col = "black",
                             min.dist = 0.3,
                             spread = 1,
                             trajectory_coords = NULL,
                             trajectory_col = NULL,
                             lwd = 2)
{
  suppressPackageStartupMessages({
    require(plotly)
    require(scales)
    require(KernSmooth)
    require(Seurat)
    require(dplyr)
  })

  # Ensure feature is in the correct assay
  DefaultAssay(seurat_object) = assay

  cat("Calculating 3D components\n")
  temp = Seurat::RunUMAP(object = seurat_object,
                         dims = dims,
                         n.components = 3,
                         min.dist = min.dist,
                         spread = spread,
                         verbose = F,
                         assay = assay,
                         reduction = reduction)

  cat("Extracting feature data\n")
  feature_data = Seurat::FetchData(temp,
                                   vars = gene)

  # Create a data frame with the UMAP coordinates and feature expression
  df = data.frame(umap1 = temp@reductions$umap@cell.embeddings[, 1],
                  umap2 = temp@reductions$umap@cell.embeddings[, 2],
                  umap3 = temp@reductions$umap@cell.embeddings[, 3],
                  feature_expression = feature_data[, gene])

  # Normalize feature expression for color scaling
  df$feature_expression = scales::rescale(df$feature_expression,
                                           to = c(0, 1))

  cat("Constructing 3D Feature Plot\n")
  p = plotly::plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3,
                      color = ~feature_expression,
                      colors = cols,
                      marker = list(size = pt.size),
                      stroke = stroke,
                      strokes = "black",
                      type = "scatter3d",
                      mode = "markers") %>%
    plotly::layout(scene = list(xaxis = list(title = "UMAP_1"),
                                yaxis = list(title = "UMAP_2"),
                                zaxis = list(title = "UMAP_3"))) %>%
    plotly::layout(plot_bgcolor = bg_col) %>%
    plotly::layout(paper_bgcolor = bg_col)

  if (!is.null(trajectory_coords)) {
    trajectory_col = trajectory_col %||% "black"
    for (i in 1:nrow(trajectory_coords)) {
      p = p %>% plotly::add_trace(x = as.vector(t(trajectory_coords[i, c("source_dim1", "target_dim1")])),
                                  y = as.vector(t(trajectory_coords[i, c("source_dim2", "target_dim2")])),
                                  z = as.vector(t(trajectory_coords[i, c("source_dim3", "target_dim3")])),
                                  line = list(color = I(trajectory_col),
                                              width = lwd),
                                  type = "scatter3d",
                                  mode = "lines",
                                  showlegend = FALSE, inherit = F)
    }
  }

  return(p)
  gc()
}
