#' Dot plot of relative gene expression
#'
#' This function wraps the Seurat DotPlot function with ggplot2 functions, for added customizability.
#' @param seurat_object A Seurat object.
#' @param genes A vector of genes to plot.
#' @param cols A vector of colors to use for low and high values, respectively.
#' @param border Whether or not to draw a border to the dots.
#' @param title A title for the plot.
#' @param portrait Whether or not to plot genes on the y-axis
#' @param font Font to use for text.
#' @param plot_missing Logical; whether genes not found in the Seurat object should still be included in the plot.
#' @export
#' @examples
#' KG_dotplot(B_cell_dataset, genes = c("CD19", "CD27", "CD38", "CR2"))
KG_dotplot = function(seurat_object, genes, cols = c("blue4", "red2"),
                      border = TRUE,
                      title = NULL,
                      portrait = TRUE,
                      font = "Arial",
                      plot_missing = FALSE) {

  suppressPackageStartupMessages({
    require(Seurat)
    require(ggplot2)
    require(stats)
    require(dplyr)
    require(S4Vectors)
    require(plotly)})

  # Extract values from Seurat's DotPlot function
  df = Seurat::DotPlot(object = seurat_object, features = genes)$data

  if(plot_missing){
    missing_genes = setdiff(genes, unique(df$features.plot))

    missing_df = S4Vectors::expand.grid(features.plot = missing_genes,
                                        id = unique(df$id)) %>%
      plotly::mutate(avg.expr.scaled = 0,
                    pct.exp = 0)

    df = dplyr::bind_rows(df, missing_df)
  }

  # Construct the plot
  plot = ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(features.plot, dplyr::desc(features.plot)), y = reorder(id, dplyr::desc(id)))) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, fill = avg.exp.scaled),
               color = if(border) "black" else "transparent",
               shape = 21) +
    ggplot2::scale_size("% detected", range = c(0, 9)) +
    ggplot2::scale_fill_gradient(low = cols[1], high = cols[2],
                        guide = ggplot2::guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Relative\nexpression") +
    ggplot2::ylab("Cluster") +
    ggplot2::xlab("") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 18, family = font),
          axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, color = "black", family = font),
          axis.text.y = ggplot2::element_text(size = 12, color = "black", family = font),
          axis.title = ggplot2::element_text(size = 14, family = font)) +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(color = "black")))

  if(portrait){
    plot = plot + ggplot2::coord_flip()
  }

  return(plot)
}
