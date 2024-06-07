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
#' @param axes Logical; whether or not to plot the axes.
#' @param heading Logical; whether or not to include a heading (with the name of the isotype plotted).
#' @param isotype Slow in Seurat object where the isotype information is stored.
#' @param combine Logical; whether or not to print all plots on a single page.
#' @param ncol Number of columns to organize the plots in when printed.
#' @export
#' @examples
#' KG_isotype_dimplot(B_cell_dataset)

KG_isotype_dimplot = function(seurat_object,
                              cells = NULL,
                              reduction = "umap",
                              dims = c(1, 2),
                              group.by = NULL,
                              label = TRUE,
                              label.size = 7,
                              pt.size = 2,
                              stroke = 0.1,
                              legend = FALSE,
                              axes = TRUE,
                              heading = TRUE,
                              isotype = "isotype",
                              combine = TRUE,
                              ncol = 2
){

  suppressPackageStartupMessages(c(
    require(Seurat),
    require(SeuratObject),
    require(ggplot2),
    require(patchwork),
    require(dplyr)))

  # Set up the basics
  reduction = reduction %||% SeuratObject::DefaultDimReduc(seurat_object)
  seurat_object[['ident']] = factor(Seurat::Idents(seurat_object))
  group.by = group.by %||% 'ident'
  cells = cells %||% colnames(seurat_object)

  # Extract embeddings
  df = data.frame(
    Dim1 = Seurat::Embeddings(seurat_object[[reduction]])[cells, dims[1]],
    Dim2 = Seurat::Embeddings(seurat_object[[reduction]])[cells, dims[2]],
    Identity = factor(seurat_object[[group.by]][cells, ]),
    Isotype = seurat_object[[isotype]][cells, ]
  )

  # Isotypes, colors, and dimension names
  isotypes = c("IgD", "IgM", "IgG", "IgA", "IgE")
  cols = c("beige", "cornflowerblue", "brown4", "tan2", "lightgreen")

  if(reduction == "pca"){
    dim_names = c("PC1", "PC2")
  } else{dim_names = c(paste(toupper(reduction), "1", sep = ""),
                       paste(toupper(reduction), "2", sep = ""))

  }
  # Create the initial plots
  plots = list()
  for(i in seq_along(isotypes)){

    ## Create a highlighting condition for each isotype
    df$highlight = as.character(df$Isotype %in% isotypes[i])
    df = df[order(df$highlight), ]

    ## Construct plots
    plots[[i]] = ggplot2::ggplot(df, aes(x = Dim1, y = Dim2)) +
      ggplot2::geom_point(shape = 21, size = pt.size,
                          stroke = stroke,
                          color = "black",
                          aes(fill = highlight)) +
      ggplot2::scale_fill_manual(values = c("TRUE" = cols[i],
                                            "FALSE" = "lightgrey"), na.value = "lightgrey") +
      ggplot2::theme_classic() +
      ggplot2::ggtitle(isotypes[i]) +
      ggplot2::xlab(dim_names[1]) +
      ggplot2::ylab(dim_names[2]) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Arial"),
                     axis.title = element_text(hjust = 0.5, size = 20, family = "Arial"),
                     axis.text = element_text(size = 14, family = "Arial"))
  }

  # Label identity classes?
  if(label){
    names(df)[1:2] = c("Dim1", "Dim2")
    group_centers = df %>%
      dplyr::group_by(Identity) %>%
      dplyr::summarize(center_x = mean(Dim1),
                       center_y = mean(Dim2),
                       .groups = 'drop')

    for(i in 1:length(plots)){
      plots[[i]] = plots[[i]] +
        ggplot2::geom_text(data = group_centers,
                           aes(x = center_x,
                               y = center_y,
                               label = Identity),
                           size = label.size)
    }
  }
  # Legend?
  if(!legend){
    for(i in 1:length(plots)){
      plots[[i]] = plots[[i]] +
        ggplot2::theme(legend.position = "none")
    }
  }

  # Axes?
  if(!axes){
    for(i in 1:length(plots)){
      plots[[i]] = plots[[i]] +
        Seurat::NoAxes()
    }
  }

  # Heading?
  if(!heading){
    for(i in 1:length(plots)){
      plots[[i]] = plots[[i]] +
        ggplot2::theme(plot.title = element_blank())
    }
  }
  if(combine){
    return(patchwork::wrap_plots(plots, ncol = ncol))
  } else {
    return(plots)
  }
}
