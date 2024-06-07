#' Plot percentage distribution of isotypes by cluster as pie charts
#'
#' This function plots the distribution of Ig isotypes for clusters in a Seurat object.
#' @param seurat_object A Seurat object
#' @param MGA_only Whether or not to only plot IgM, IgG, and IgA isotypes
#' @param cluster Name of cluster(s) for which to generate pie charts
#' @param combine Whether or not plots should be combined into a single or separate plots
#' @param ncol The number of columns for the layout of the plot
#' @param text Logical; whether or not to print the percentage values for each slice.
#' @export
#' @examples
#' KG_isotype_pie(B_cell_dataset, cluster = "0")
KG_isotype_pie = function(seurat_object,
                          MGA_only = TRUE,
                          cluster = NULL,
                          combine = TRUE,
                          ncol = 2,
                          text = TRUE) {

  suppressPackageStartupMessages({
    require(gridExtra)
    require(ggplot2)})

  if(is.null(cluster)) {
    stop("No cluster specified")
  } else {
    # Calculate percentages for each cluster
    percent_df = data.frame() # Empty data frame to store percentages in

    if(MGA_only == TRUE){
      tab = table(Seurat::Idents(seurat_object),
                  seurat_object$isotype)[, c(1, 3, 4)]

      for(i in seq_len(nrow(tab))) {
        percentages = c(
          IgM = sum(tab[i, 1]) / sum(tab[i, ]) * 100,
          IgG = sum(tab[i, 2]) / sum(tab[i, ]) * 100,
          IgA = sum(tab[i, 3]) / sum(tab[i, ]) * 100
        )

        # Add to results data frame
        percent_df = rbind(percent_df,
                           cbind(Cluster = rownames(tab)[i],
                                 round(t(as.data.frame(percentages)),
                                       1)))
      }

      percent_df = percent_df[percent_df$Cluster %in% c(cluster),]
      plot_list = list() # Initialize an empty list to store plots

      # Create a pie chart for each cluster
      for(i in seq_len(nrow(percent_df))) {
        df_for_plot = data.frame(
          Isotype = factor(c("IgM", "IgG", "IgA"),
                           levels = rev(c("IgM", "IgG", "IgA"))),
          Value = as.integer(c(percent_df$IgM[i], percent_df$IgG[i],
                               percent_df$IgA[i])),
          Labels = paste(c(percent_df$IgM[i], percent_df$IgG[i],
                           percent_df$IgA[i]), "%", sep = "")
        )

        # Calculate the cumulative percentages for label placement
        df_for_plot$Position = cumsum(as.integer(df_for_plot$Value)) - 0.5 * as.integer(df_for_plot$Value)

        p = ggplot2::ggplot(df_for_plot, ggplot2::aes(x = "", y = Value, fill = Isotype)) +
          ggplot2::geom_bar(width = 1, stat = "identity", colour = "black", stroke = 0.5) +
          ggplot2::coord_polar("y", start = 0) +
          ggplot2::scale_fill_manual(values = c("IgM" = "cornflowerblue",
                                                "IgG" = "brown4",
                                                "IgA" = "tan1"),
                                     breaks = c("IgM", "IgG", "IgA")) +
          ggplot2::geom_text(ggplot2::aes(y = Position,
                                 label = Labels),
                             color = "white", size = 5) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste("Cluster", percent_df$Cluster[i])) +
          ggplot2::theme(legend.position = "right",
                         plot.title = ggplot2::element_text(hjust = 0.5))

        plot_list[[i]] = p
      }
    } else{
      tab = table(Seurat::Idents(seurat_object),
                  seurat_object$isotype)

      for(i in seq_len(nrow(tab))) {
        percentages = c(
          IgM = sum(tab[i, 1]) / sum(tab[i, ]) * 100,
          IgD = sum(tab[i, 2]) / sum(tab[i, ]) * 100,
          IgG = sum(tab[i, 3]) / sum(tab[i, ]) * 100,
          IgA = sum(tab[i, 4]) / sum(tab[i, ]) * 100,
          IgE = sum(tab[i, 5]) / sum(tab[i, ]) * 100
        )

        # Add to results data frame
        percent_df = rbind(percent_df,
                           cbind(Cluster = rownames(tab)[i],
                                 round(t(as.data.frame(percentages)), 1)))
      }

      percent_df = percent_df[percent_df$Cluster %in% c(cluster),]
      plot_list = list() # Initialize an empty list to store plots

      # Create a pie chart for each cluster
      for(i in seq_len(nrow(percent_df))) {
        df_for_plot = data.frame(
          Isotype = factor(c("IgM", "IgD", "IgG", "IgA", "IgE"),
                           levels = rev(c("IgM", "IgD", "IgG", "IgA", "IgE"))),
          Value = as.integer(c(percent_df$IgM[i], percent_df$IgD[i],
                               percent_df$IgG[i], percent_df$IgA[i],
                               percent_df$IgE[i])),
          Labels = paste(c(percent_df$IgM[i], percent_df$IgD[i],
                           percent_df$IgG[i], percent_df$IgA[i],
                           percent_df$IgE[i]),
                         "%", sep = "")
        )

        # Calculate the cumulative percentages for label placement
        df_for_plot$Position = cumsum(as.integer(df_for_plot$Value)) - 0.5 * as.integer(df_for_plot$Value)

        p = ggplot2::ggplot(df_for_plot, ggplot2::aes(x = "", y = Value, fill = Isotype)) +
          ggplot2::geom_bar(width = 1, stat = "identity", colour = "black") +
          ggplot2::coord_polar("y", start = 0) +
          ggplot2::scale_fill_manual(values = c("IgM" = "cornflowerblue",
                                                "IgD" = "wheat2",
                                                "IgG" = "brown4",
                                                "IgA" = "tan1",
                                                "IgE" = "lightgreen"),
                                     breaks = c("IgM", "IgD", "IgG", "IgA", "IgE")) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste("Cluster", percent_df$Cluster[i])) +
          ggplot2::theme(legend.position = "right",
                         plot.title = ggplot2::element_text(hjust = 0.5))
        if(text){
          p = p +
            ggplot2::geom_text(ggplot2::aes(y = Position,
                                   label = Labels),
                               color = "white", size = 5)
        }
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
