#' Make a nice feature plot!
#'
#' This function constructs nice, clean feature plots from dimensionally reduced data points in a Seurat object.
#' @param seurat_object A Seurat object.
#' @param genes A gene, or vector of genes, for which expression should be plotted.
#' @param reduction The dimensionality reduction to use.
#' @param dims Dimensions to plot. By default, plots the first two dimensions.
#' @param label Whether to plot group identifier labels.
#' @param label.size Size of the group identifier labels.
#' @param pt.size Size of data points.
#' @param stroke Size of border to data points. NOTE: when exportet as PDF, borders can alter the appearance of the plot. To avoid this, export externally as PNG file.
#' @param legend Whether or not to plot the legend/key to the group identifier.
#' @param cols A vector of colors to use for the shading of data points. Should be a vector of three colors, for low, mid, and high expression.
#' @param axes Logical; whether axes should be shown.
#' @param legend Logical; whether to plot the legend.
#' @export
#' @examples
#' KG_featureplot(B_cell_dataset, "CD27")
KG_featureplot = function(seurat_object,
                          genes,
                          reduction = "umap",
                          dims = c(1, 2),
                          pt.size = 2,
                          cols = c("#FEFDE2", "#FC8C40", "#950028"),
                          stroke = 0.05,
                          axes = TRUE,
                          legend = TRUE){
  suppressPackageStartupMessages({
    require(Seurat)
    require(dplyr)
    require(gridExtra)
    require(plotly)
    require(ggplot2)})

  plot_list = list()

  for (i in seq_along(genes)){
    message(genes[i])
    plot_list[[i]] = local({
      i = i

      a = data.frame(dim1 = Seurat::Embeddings(seurat_object[[reduction]])[, dims[1]],
                     dim2 = Seurat::Embeddings(seurat_object[[reduction]])[, dims[2]],
                     Seurat::FetchData(seurat_object, vars = genes[i]))

      axis_labs = c(paste0(toupper(reduction), "_", as.character(dims[1])),
                    paste0(toupper(reduction), "_", as.character(dims[2])))

      names(a)[3] = "Gene"
      a = plotly::arrange(a, desc(Gene))
      p1 = ggplot2::ggplot(data = a %>% plotly::arrange(Gene),
                           aes(x = dim1, y = dim2)) +
        ggplot2::geom_point(aes(fill = Gene),
                            shape = 21,
                            size = pt.size,
                            stroke = stroke) +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_gradient2(low = cols[1],
                                      mid = cols[2],
                                      high = cols[3],
                                      midpoint = 1,
                                      na.value = "lightgrey",
                                      limits = c(0.5, max(a$Gene))) +
        ggplot2::labs(fill = "Expression") +
        ggplot2::guides(fill = ggplot2::guide_colorbar(ticks.colour = NA)) +
        ggplot2::ggtitle(genes[i]) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 18),
                       legend.title = ggplot2::element_text(size = 14),
                       legend.text = ggplot2::element_text(size = 11),
                       axis.title = ggplot2::element_text(size = 14),
                       axis.text = ggplot2::element_text(size = 11)) +
        ggplot2::xlab(axis_labs[1]) +
        ggplot2::ylab(axis_labs[2])

      if(axes == FALSE){
        p1 = p1 + Seurat::NoAxes()
      }

      if(legend == FALSE){
        p1 = p1 + ggplot2::theme(legend.position = "none")
      }

      return(p1)
    })
  }

  ncol_plot = ifelse(length(genes) > 1, 2, 1)
  do.call(gridExtra::grid.arrange, c(plot_list, ncol = ncol_plot))
}
