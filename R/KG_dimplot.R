#' Make a nice dimplot
#'
#' This function constructs a nice, clean dimplot from dimensionally reduced data points in a Seurat object.
#' @param seurat_object A Seurat object.
#' @param cells List of cells to plot. If NULL, all cells are plotted.
#' @param reduction The dimensionality reduction to use.
#' @param dims Dimensions to plot. By default, plots the first two dimensions.
#' @param group.by Group identifier to use for shading.
#' @param label Whether to plot group identifier labels.
#' @param label.size Size of the group identifier labels.
#' @param pt.size Size of data points.
#' @param stroke Size of border to data points. NOTE: when exportet as PDF, borders can alter the appearance of the plot. To avoid this, export externally as PNG file.
#' @param legend Whether or not to plot the legend/key to the group identifier.
#' @param cols A vector of colors to use for the shading of data points.
#' @param trajectory_coords A vector containing coordinates for drawing trajectory lines.
#' @param trajectory_col What color to use for the trajectory.
#' @param pub_ready Layout option. When TRUE, produces a more clean plot, with shorter axis lines and no ticks.
#' @export
#' @examples
#' KG_dimplot(B_cell_dataset)
KG_dimplot = function(seurat_object,
                      cells = NULL,
                      reduction = NULL,
                      dims = c(1, 2),
                      group.by = NULL,
                      label = FALSE,
                      label.size = 4,
                      pt.size = 3,
                      stroke = 0.2,
                      legend = TRUE,
                      cols = NULL,
                      trajectory_coords = NULL,
                      trajectory_col = NULL,
                      lwd = 0.5,
                      pub_ready = FALSE) {
  suppressPackageStartupMessages(c(require(Seurat),
                                   require(SeuratObject),
                                   require(ggplot2),
                                   require(viridis),
                                   require(scales),
                                   require(colorspace),
                                   require(grDevices),
                                   require(grid)))

  # Define reductions, identities, and cells
  reduction = reduction %||% SeuratObject::DefaultDimReduc(seurat_object)
  seurat_object[['ident']] = factor(Seurat::Idents(seurat_object))
  group.by = group.by %||% 'ident'
  cells = cells %||% colnames(seurat_object)

  # Extract cell embeddings
  a = data.frame(Seurat::Embeddings(seurat_object[[reduction]])[cells, dims[1]],
                 Seurat::Embeddings(seurat_object[[reduction]])[cells, dims[2]],
                 seurat_object[[group.by]][cells, ])
  a[, 3] = factor(a[, 3])
  dims = paste0(Key(seurat_object[[reduction]]), dims)

  names(a) = c(dims[1],
               dims[2],
               group.by)

  # Colors
  if(is.null(cols)){
    cols = scales::hue_pal()(length(levels(seurat_object[[group.by, drop = T]])))
  }
  names(a)[3] = "groups"
  p = ggplot2::ggplot(a, aes(x = a[,1], y = a[,2], fill = groups)) +
    ggplot2::geom_point(shape = 21, size = pt.size, stroke = stroke) +
    ggplot2::theme_classic() +
    #labs(fill = legend_title) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = names(a)[1], y = names(a)[2])

  if(!is.null(trajectory_coords)){
    trajectory_col = trajectory_col %||% "black"
    p = p +
      ggplot2::geom_segment(data = trajectory_coords, aes(x = source_dim1, xend = target_dim1,
                                                          y = source_dim2, yend = target_dim2),
                            color = trajectory_col, size = lwd,
                            inherit.aes = FALSE)
  }

  if(pub_ready){
    p = p +
      # Set new axis lines, with arrows
      ggplot2::geom_segment(aes(x = min(a[,1]) - 0.25 , y = min(a[,2]) - 0.25,
                                xend = min(a[,1]) + 2.5, yend = min(a[,2]) - 0.25),
                            arrow = grid::arrow(length = unit(.3, 'cm'), type = "closed")) +
      ggplot2::geom_segment(aes(x = min(a[,1]) - 0.25, y = min(a[,2]) - 0.25,
                                xend = min(a[,1]) - 0.25, yend = min(a[,2]) + 2.5),
                            arrow = grid::arrow(length = unit(.3, 'cm'), type = "closed")) +

      # Update theme to set text size and remove the original axis lines
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 11),
                     axis.title.x = element_text(hjust = 0.06, vjust = 5, size = 14),
                     axis.title.y = element_text(hjust = 0.06, vjust = -4, angle = 90, size = 14),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.line = element_blank(),
                     # Since manually setting axis lines with geom_segment changes the margins of the plot,
                     # reset the margins to visually equal to when pub_ready = FALSE
                     plot.margin = margin(12, 12, 12, 12))
  } else{
    p = p + ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18),
                           legend.title = element_text(size = 14),
                           legend.text = element_text(size = 11),
                           axis.title.x = element_text(hjust = 0.5, size = 14),
                           axis.title.y = element_text(hjust = 0.5, angle = 90, size = 14),
                           axis.text = element_text(size = 10))
  }

  # Legend
  if(!legend) {
    p = p + ggplot2::theme(legend.position = "none")
  }

  # Labels
  if(label == TRUE){
    b = a
    names(b)[1:2] = c("Dim1", "Dim2")

    ## Calculate cluster centroids
    group_centers = b %>%
      dplyr::group_by(groups) %>%
      dplyr::summarize(center_x = mean(Dim1), center_y = mean(Dim2))

    ## Add labels to plot
    p = p + ggplot2::geom_text(data = group_centers,
                               aes(x = center_x,
                                   y = center_y,
                                   label = groups,
                                   family = "Arial"),
                               size = label.size)
  }

  return(p)
}
