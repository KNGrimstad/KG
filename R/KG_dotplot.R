#' Dot plot of relative gene expression
#'
#' This function wraps the Seurat DotPlot function with ggplot2 functions, for added customizability.
#' @param seurat_object A Seurat object.
#' @param genes A vector of genes to plot.
#' @param cols A vector of colors to use for low and high values, respectively.
#' @param border Whether or not to draw a border to the dots.
#' @param title A title for the plot.
#' @param portrait Whether or not to plot genes on the y-axis
#' @export
#' @examples
#' KG_dotplot(B_cell_dataset, genes = c("CD19", "CD27", "CD38", "CR2"))
KG_dotplot = function(seurat_object, genes, cols = c("blue4", "red2"),
                      border = TRUE, title = NULL, portrait = TRUE, font = "arial") {
  suppressPackageStartupMessages(c(require(Seurat),
                                   require(ggplot2)))

  # Extract values from Seurat's DotPlot function
  df = Seurat::DotPlot(object = seurat_object, features = genes)$data

  # Construct the plot
  plot = ggplot2::ggplot(df, aes(x = reorder(features.plot, desc(features.plot)), y = reorder(id, desc(id)))) +
    ggplot2::geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
               color = if(border) "black" else "transparent",
               shape = 21) +
    ggplot2::scale_size("% detected", range = c(0, 9)) +
    ggplot2::scale_fill_gradient(low = cols[1], high = cols[2],
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Relative\nexpression") +
    ggplot2::ylab("Cluster") + xlab("") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 18, family = font),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black", family = font),
          axis.text.y = element_text(size = 12, color = "black", family = font),
          axis.title = element_text(size = 14, family = font)) +
    ggplot2::guides(size = guide_legend(override.aes = list(color = "black")))

  if(portrait){
    plot = plot + ggplot2::coord_flip()
  }

  return(plot)
}
