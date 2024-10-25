#' 3D UMAP for Seurat object
#'
#' This function computes and plots cells in three dimensions.
#' @param seurat_object A Seurat object.
#' @param reduction Name of the dimensional reduction method to use as input, "pca" by default.
#' @param dims The range of dimensions to use from linear dimensional reduction method.
#' @param assay Name of the Seurat object assay to use.
#' @param pt.size Controls the size of plotted points/cells.
#' @param stroke Size of border to data points.
#' @param cols Determines the colors to use for cells. If NULL, uses default colors.
#' @param bg.col Controls the background color of the plot.
#' @param min.dist Controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimize more accurately the local structure.
#' @param spread Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.
#' @param group.by Controls how cells should be grouped/colored. By default, uses the default identity of cells in the Seurat object.
#' @param trajectory_coords A vector containing coordinates for drawing trajectory lines.
#' @param trajectory_col What color to use for the trajectory.
#' @param lwd Line width for trajectory.
#' @export
#' @examples
#' KG_3DPlot(B_cell_dataset)
KG_3DPlot = function(seurat_object,
                     dims = 1:35,
                     reduction = "pca",
                     assay = "RNA",
                     pt.size = 7,
                     stroke = 1,
                     cols = NULL,
                     bg.col = "black",
                     min.dist = 0.3,
                     spread = 1,
                     group.by = NULL,
                     trajectory.coords = NULL,
                     trajectory.col = NULL,
                     lwd = 2){

  suppressPackageStartupMessages({
    require(plotly)
    require(scales)
    require(KernSmooth)
    require(Seurat)
    require(dplyr)})

  # Calculate UMAP for three components
  if (!is.null(group.by)){
    Seurat::Idents(seurat_object) = seurat_object[[group.by, drop = TRUE]]
  }
  cat("Calculating 3D components\n")
  temp = Seurat::RunUMAP(object = seurat_object,
                         dims = dims,
                         n.components = 3,
                         min.dist = min.dist,
                         spread = spread,
                         verbose = F,
                         assay = assay,
                         reduction = reduction)

  # Extract data
  df = data.frame(umap1 = temp@reductions$umap@cell.embeddings[, 1],
                  umap2 = temp@reductions$umap@cell.embeddings[, 2],
                  umap3 = temp@reductions$umap@cell.embeddings[, 3],
                  cell = Seurat::Idents(temp))

  # Select color palette
  if(is.null(cols)){
    cols = scales::hue_pal()(length(unique(df$cell)))
  }

  # Plot in 3D
  cat("Constructing 3D view\n")
  p = plotly::plot_ly(df,
                      x = ~umap1,
                      y = ~umap2,
                      z = ~umap3,
                      color = ~cell,
                      marker = list(size = pt.size),
                      colors = cols,
                      stroke = stroke,
                      strokes = "black") %>%

    plotly::add_markers() %>%
    plotly::layout(scene = list(xaxis = list(title = "UMAP_1"),
                                yaxis = list(title = "UMAP_2"),
                                zaxis = list(title = "UMAP_3"))) %>%
    plotly::layout(plot_bgcolor = bg.col) %>%
    plotly::layout(paper_bgcolor = bg.col)

  # Add trajectory
  if(!is.null(trajectory.coords)){
    trajectory.col = trajectory.col %||% "black"

    for(i in 1:nrow(trajectory.coords)){
      p = p %>%
        plotly::add_trace(
          x = as.vector(t(trajectorycoords[i, c("source_dim1", "target_dim1")])),
          y = as.vector(t(trajectory.coords[i, c("source_dim2", "target_dim2")])),
          z = as.vector(t(trajectory.coords[i, c("source_dim3", "target_dim3")])),
          #color = trajectory_col,
          line = list(color = I(trajectory.col),
                      width = lwd),
          type = "scatter3d",
          mode = "lines",
          showlegend = FALSE,
          inherit = F)
    }
  }
  return(p)
  gc()
}
