#' Plot percentage distribution of isotype subtypes by cluster as pie charts
#'
#' This function plots the distribution of Ig subtypes of IgG (default) or IgA for clusters in a Seurat object.
#' @param seurat_object A Seurat object
#' @param subtype Whether to plot IgG or IgA subtypes
#' @param c_call The name of the meta data slot in which the heavy chain C-gene calls are stored
#' @param cluster The name(s) of then cluster(s) for which to generate pie charts
#' @param combine Whether or not plots should be combined into a single or separate plots
#' @param ncol The number of columns for the layout of the plot
#' @export
#' @examples
#' KG_subtype_pie(B_cell_dataset, c_call = "c_call", cluster = "0")
KG_subtype_pie = function(seurat_object,
                           subtype = "IgG",
                           c_call = "c_call",
                           cluster = NULL,
                           combine = TRUE,
                           ncol = 2) {

  suppressPackageStartupMessages({
    require(gridExtra)
    require(ggplot2)})

  if(is.null(cluster)) {
    stop("No cluster specified")
  } else {
    # Check that the specified column (c_call) exists in the Seurat object
    if (!(c_call %in% names(seurat_object@meta.data))) {
      stop(paste("Meta data slot", c_call, "not found in the Seurat object"))
    }
    # Calculate percentages for each cluster
    percent_df <- data.frame() # Empty data frame to store percentages in

    tab = table(Seurat::Idents(seurat_object),
                seurat_object[[c_call, drop = T]])

    if(subtype == "IgG"){
      tab = tab[, c(5, 6, 7, 8)]
    } else if(subtype == "IgA"){
      tab = tab[, c(1, 2)]
    } else{
      stop("Invalid specification of isotype: must be IgG or IgA")
    }

    for(i in seq_len(nrow(tab))) {
      if((subtype == "IgG")){
        percentages = c(
          IgG1 = sum(tab[i, 1]) / sum(tab[i, ]) * 100,
          IgG2 = sum(tab[i, 2]) / sum(tab[i, ]) * 100,
          IgG3 = sum(tab[i, 3]) / sum(tab[i, ]) * 100,
          IgG4 = sum(tab[i, 4]) / sum(tab[i, ]) * 100
        )
      } else{
        percentages = c(
          IgA1 = sum(tab[i, 1]) / sum(tab[i, ]) * 100,
          IgA2 = sum(tab[i, 2]) / sum(tab[i, ]) * 100
        )
      }

      # Add to results data frame
      percent_df = rbind(percent_df, cbind(Cluster = rownames(tab)[i],
                                           round(t(as.data.frame(percentages)), 1)))
    }

    percent_df = percent_df[percent_df$Cluster %in% c(cluster),]
    plot_list = list() # Initialize an empty list to store plots

    # Create a pie chart for each cluster
    for(i in seq_len(nrow(percent_df))) {
      if(subtype == "IgG"){
        df_for_plot = data.frame(
          "Subtype" = factor(c("IgG1", "IgG2", "IgG3", "IgG4"),
                           levels = c("IgG4", "IgG3", "IgG2", "IgG1")),
          "Value" = as.integer(c(percent_df$IgG1[i], percent_df$IgG2[i],
                               percent_df$IgG3[i], percent_df$IgG4[i])),
          "Labels" = paste(c(percent_df$IgG1[i], percent_df$IgG2[i],
                           percent_df$IgG3[i], percent_df$IgG4[i]),
                         "%", sep = "")
        )

      } else {
        df_for_plot = data.frame(
          "Subtype" = factor(c("IgA1", "IgA2"),
                           levels = c("IgA2", "IgA1")),
          "Value" = as.integer(c(percent_df$IgA1[i], percent_df$IgA2[i])),
          "Labels" = paste(c(percent_df$IgA1[i], percent_df$IgA2[i])),
          "%", sep = "")
      }

      # Calculate the cumulative percentages for label placement
      df_for_plot$Position = cumsum(as.integer(df_for_plot$Value)) - 0.5 * as.integer(df_for_plot$Value)

      if(subtype == "IgG"){
        p = ggplot2::ggplot(df_for_plot, ggplot2::aes(x = "", y = Value, fill = Subtype)) +
          ggplot2::geom_bar(width = 1, stat = "identity", colour = "black") +
          ggplot2::coord_polar("y", start = 0) +
          ggplot2::scale_fill_manual(values = c("IgG1" = "brown4",
                                                "IgG2" = "indianred3",
                                                "IgG3" = "indianred1",
                                                "IgG4" = "darksalmon"),
                                     breaks = c("IgG1", "IgG2", "IgG3", "IgG4")) +
          ggplot2::geom_text(aes(y = Position, label = Labels), color = "white", size = 5) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste("Cluster", percent_df$Cluster[i])) +
          ggplot2::theme(legend.position = "right",
                         plot.title = ggplot2::element_text(hjust = 0.5))

        plot_list[[i]] = p
      } else {
        p = ggplot2::ggplot(df_for_plot, ggplot2::aes(x = "", y = Value, fill = Subtype)) +
          ggplot2::geom_bar(width = 1, stat = "identity", colour = "black") +
          ggplot2::coord_polar("y", start = 0) +
          ggplot2::scale_fill_manual(values = c("IgA1" = "tan1",
                                                "IgA2" = "tan3"),
                                     breaks = c("IgA1", "IgA2")) +
          ggplot2::geom_text(aes(y = Position, label = Labels), color = "white", size = 5) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste("Cluster", percent_df$Cluster[i])) +
          ggplot2::theme(legend.position = "right",
                         plot.title = ggplot2::element_text(hjust = 0.5))

        plot_list[[i]] = p
      }
    }

    # Combine plots if required
    if(combine && length(plot_list) > 1) {
      gridExtra::grid.arrange(grobs = plot_list,
                              ncol = ncol)
    } else {
      return(plot_list)
    }
  }
}
