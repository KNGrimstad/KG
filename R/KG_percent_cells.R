#' Plot percentage of cells for each cluster.
#'
#' This function plots the percentage of cells, out of the entire dataset, that each cluster contains.
#' @param seurat_object A Seurat object
#' @param group.by The identity column to use from the Seurat object.
#' @param plot How to plot the results. Options include "table", "pie", "bar".
#' @param cols Vector of colors to use.
#' @param label Logical; whether or not to add the percentages as text albels in the plot. Irrelevant for tables.
#' @param label_size The size of labels.
#' @param label_col Color for the labels.
#' @param legend Logical; whether or not to include legend.
#' @export
#' @examples
#'KG_percent_cells(B_cell_dataset)
KG_percent_cells = function(seurat_object,
                            group.by = NULL,
                            plot = "pie",
                            cols = NULL,
                            label = T,
                            label_size = 3,
                            label_col = "black",
                            legend = T){
  suppressPackageStartupMessages(c(require(gridExtra),
                                   require(ggplot2),
                                   require(scales),
                                   require(Seurat)))
  # Define stuff
  seurat_object[['ident']] = factor(Seurat::Idents(seurat_object))
  group.by = group.by %||% 'ident'
  cols = cols %||% scales::hue_pal()(length(levels(seurat_object[[group.by, drop = T]])))

  df = data.frame(round(prop.table(table(seurat_object[[group.by, drop = T]])) * 100, 1))
  names(df) = c("Cluster", "Percentage")

  df_for_plot = data.frame(Cluster = factor(df$Cluster, levels = levels(df$Cluster)),
                           Percentage = c(round(df$Percentage, 1)))

  if(plot == "pie"){
    labels = vector()
    for(i in 1:length(df$Cluster)){
      labels[i] = paste(round(df$Percentage[i], 1), "%", sep = "")
    }
    df_for_plot$Labels = labels
    df_for_plot$Position = cumsum(df_for_plot$Percentage) - 0.25 * df_for_plot$Percentage

    p = ggplot2::ggplot(df_for_plot, aes(x = "", y = rev(Percentage), fill = Cluster)) +
      ggplot2::geom_bar(width = 1, stat = "identity", colour = "black") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_manual(values = rev(cols)) +
      ggplot2::theme_void()

    if(legend){
      p = p +
        ggplot2::theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5))
    } else {
      p = p +
        ggplot2::theme(legend.position = "none",
                       plot.title = element_text(hjust = 0.5))
    }
    if(label){
      p = p +
        ggplot2::geom_text(aes(x = 1, label = rev(Labels)),
                           color = label_col,
                           size = label_size,
                           position = ggplot2::position_stack(vjust = 0.5))
    }

  } else if(plot == "donut"){
    labels = vector()
    for(i in 1:nrow(df)){
      labels[i] = paste(round(df$Percentage[i], 1), "%", sep = "")
    }
    df_for_plot$Labels = labels # NOTE: labels are currently not supported
    df_for_plot$Position = cumsum(df_for_plot$Percentage) - 0.5 * df_for_plot$Percentage


    p = ggplot2::ggplot(df_for_plot, aes(x = "", y = rev(Percentage), fill = Cluster)) +
      ggplot2::geom_bar(width = 1, stat = "identity", colour = "black") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_manual(values = rev(cols)) +
      ggplot2::geom_text(aes(label = rev(Percentage),
                             x = 1.32),
                         position = ggplot2::position_stack(vjust = 0.5),
                         size = label_size,
                         color = label_col) +
      ggplot2::theme_void() +
      ggplot2::annotate("rect", xmin = -Inf, xmax = 1.105,
                        ymin = 0, ymax = Inf, fill = "black", color = "black") +
      ggplot2::annotate("rect", xmin = -Inf, xmax = 1.1,
                        ymin = 0, ymax = Inf, fill = "white") +
      ggplot2::annotate("text", x = -Inf, y = -Inf,
                        label = paste("n=",
                                      sum(table(seurat_object[[group.by, drop = T]])),
                                      sep = ""), size = 12)
    if(legend){
      p = p +
        ggplot2::theme(legend.position = "right")
    } else{
      p = p +
        ggplot2::theme(legend.position = "none")
    }


  } else if(plot == "table"){
    stop("Tables are currently not supported. \nThis feature will be added in the future.")

  } else if(plot == "bar"){
    df_for_plot$Position = df_for_plot$Percentage +0.5
    p = ggplot2::ggplot(df_for_plot, aes(x = Cluster, y = Percentage, fill = Cluster)) +
      ggplot2::geom_bar(stat = "identity", colour = "black") +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title = element_text(size = 16),
                     axis.text = element_text(size = 12),
                     axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

    if(legend){
      p = p +
        ggplot2::theme(legend.position = "right",
                       plot.title = element_text(hjust = 0.5))
    } else {
      p = p +
        ggplot2::theme(legend.position = "none",
                       plot.title = element_text(hjust = 0.5))
    }
    if(label){
      p = p +
        ggplot2::geom_text(data = df_for_plot, aes(label = paste0(Percentage, "%", sep = ""),
                                                   y = Position),
                           position = ggplot2::position_stack(),
                           size = label_size, color = label_col)
    }

  }
  return(p)
  gc()
}
