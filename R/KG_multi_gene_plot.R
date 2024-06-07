#' View feature or violin plots one by one.
#'
#' This function generates and plots individual feature or violin plots from the Seurat package for each input gene one by one. This is handy when you have a long list of genes and you want to examine the actual expression (as opposed to relative expression in e.g., heatmaps) for each gene, but cannot be bothered to print everything to a PDF file or wrap the plots.
#' @param seurat_object A Seurat object
#' @param genes A vector of genes to use for plotting.
#' @param plot Which kind of plot you want ("feature" or "violin"). Default is "feature".
#' @param pt.size Size of dots in feature plots.
#' @param feature_cols A vector of colors to use for the feature plot. Default are "lightgrey" and "red".
#' @param feature_label Logical; whether or not to plot identity labels on the feature plot.
#' @param feature_label_size Size fo the identity labels in the feature plot.
#' @param vln_cols A vector of colors to use for identities in the violin plot.
#' @param group.by Name of the identity factor in the Seurat object to use. If NULL, the current identity factor is used.
#' @param verbose Logical; text or no text?
#' @export
#' @examples
#' KG_multi_gene_plot(B_cell_dataset, genes = c("CD19", "MS4A1", "CD27", "CR2"))
KG_multi_gene_plot = function(seurat_object,
                              genes = NULL,
                              plot = "feature",
                              pt.size = 2,
                              feature_cols = c("lightgrey", "red"),
                              feature_label = FALSE,
                              feature_label_size,
                              vln_cols = NULL,
                              group.by = NULL,
                              verbose = TRUE){

  suppressPackageStartupMessages({
    require(Seurat)
    require(scales)})

  if(is.null(genes)){
    stop("You forgot to include the list of genes.")
  }

  # Identity
  if(is.null(group.by)){
    seurat_object[["ident"]] = Seurat::Idents(seurat_object)
    group.by = "ident"
  }
  Seurat::Idents(seurat_object) = seurat_object[[group.by]]

  # Violin plot colors for identities
  if(plot == "violin"){
    if(is.null(vln_cols)){
      vln_cols = scales::hue_pal()(length(unique(seurat_object[[group.by, drop = T]])))
    }
  }

  # Start loop
  for(i in 1:length(genes)){
    tryCatch({

      # Clear plotting window
      plot.new()

      # Generate plot
      if(verbose){
        cat("Showing plot ", i, " of ", length(genes),
            "\nPress [Enter] to continue, or type 'x' to exit.\n", sep = "")
      }

      # Feature plot
      if(plot == "feature"){
        plot(Seurat::FeaturePlot(seurat_object,
                                 genes[i],
                                 cols = cols,
                                 order = T,
                                 pt.size = pt.size, label = feature_label,
                                 label.size = feature_label_size) +
               Seurat::NoAxes())
      }

      # Violin plot
      if(plot == "violin"){
        plot(Seurat::VlnPlot(seurat_object,
                             genes[i],
                             cols = vln_cols))
      }

      # Pause execution
      flush.console()
      input = readline()
      if(tolower(input) == 'x'){
        cat("Exiting loop...\n")
        break
      }
    },
    error = function(e) {
      message(sprintf("An error occurred in plot %d: %s\n", i, e$message))
      if(verbose){
        cat("Press [Enter] to continue, or type 'x' to exit.\n")
        input = readline()
        if(tolower(input) == 'x') {
          cat("Exiting loop...\n")
          break
        }
      }
    })
  }
  cat("Done.\n")
}
