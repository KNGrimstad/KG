#' Bulk heatmap
#'
#' This functions draws a heatmap displaying the relative cluster-average expression of specific genes.
#' @param seurat_object A Seurat object.
#' @param genes A vector of genes for which to plot the expression.
#' @param assay The name of the assay to extract expression values from.
#' @param group_by Group variable, if NULL, uses the active identity.
#' @param show_genes Whether or not to show the gene names on the y-axis.
#' @param cols A vector of colors used for plotting the relative expression.
#' @param heading Whether or not to include a heaing to the heatmap.
#' @param title The title, if heading = TRUE.
#' @param group_bar_fill_cols A vector of colors to use for the group bars. If NULL, default colors are used.
#' @param border_col The color for the border of the tiles.
#' @param cluster_rows Whether or not genes should be clustered.
#' @param cluster_columns Whether or not groups should be clustered.
#' @param split_row_at A value, or vector of values, detailing at which row numbers a gap should be included. If NULL, no gaps.
#' @param split_column_at A value, or vector of values, detailing at which column numbers a gap should be included. If NULL, no gaps.
#' @param raster Whether or not the heatmap should be rasterized.
#' @export
#' @examples
#' # KG_heatmap(B_cell_dataset, genes = c("CD19", "CD27", "HLA-DRB1"))
KG_heatmap = function(seurat_object,
                      genes,
                      assay = "RNA",
                      group_by = NULL,
                      show_genes = TRUE,
                      cols = c("#FF00FF", "black", "#FFFF00"),
                      heading = TRUE,
                      title = NA,
                      group_bar_fill_cols = NULL,
                      border_col = "black",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      split_row_at = NULL,
                      split_column_at = NULL,
                      raster = FALSE) {
  require(ComplexHeatmap)
  require(scales)
  require(circlize)

  # Define the color function
  col_fun = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = cols)

  # Set the group_by identifier as the active identity class
    if(!is.null(group_by)){
    Idents(seurat_object) = factor(seurat_object[[group_by, drop = TRUE]])
  }

  # Calculate average expression of genes for each cluster
  avg_matrix = t(scale(t(as.matrix(AverageExpression(seurat_object, features = genes)[[assay]]))))

  # Annotations for clusters
  annotations_df = data.frame(factor(levels(Idents(seurat_object))))
  names(annotations_df) = "Cluster"
  rownames(annotations_df) = levels(Idents(seurat_object))

  # Assign colors to the clusters if not specified
  if (is.null(group_bar_fill_cols)) {
    group_bar_fill_cols = scales::hue_pal()(length(unique(annotations_df$Cluster)))
  }
  names(group_bar_fill_cols) = unique(annotations_df$Cluster)

  # Prepare annotations for ComplexHeatmap
  ha = complexHeatmap::HeatmapAnnotation(df = data.frame(Cluster = annotations_df$Cluster),
                         col = list(Cluster = group_bar_fill_cols),
                         which = "col",
                         show_annotation_name = FALSE,
                         border = TRUE,
                         gp = gpar(col = "black"),
                         show_legend = FALSE)

  # Dynamic font size
  dynamic_font = max(5, min(15, 20 - length(genes) / 5))

  # Create the heatmap
  ht = complexHeatmap::Heatmap(avg_matrix,

               # Overall layout
               border = TRUE,
               border_gp = gpar(col = "black", lty = 1, lwd = 1),
               rect_gp = gpar(col = "black", lwd = 1),
               col = col_fun,
               use_raster = raster,

               # Rows
               cluster_rows = cluster_rows,
               show_row_names = show_genes,
               row_title = NULL, # Since ComplexHeatmap manages row and column titles differently
               row_dend_reorder = TRUE,
               row_split = split_row_at,
               row_names_side = "left",
               row_names_gp = grid::gpar(fontsize = dynamic_font),

               # Columns
               cluster_columns = cluster_columns,
               column_title_gp = grid::gpar(fontsize = 15),
               column_names_gp = grid::gpar(fontsize = 15),
               show_column_names = TRUE,
               column_title = if (!is.na(title) && heading) title else NULL,
               column_dend_reorder = TRUE,
               column_split = split_column_at,

               # Legend
               name = "Average\nExpression",
               heatmap_legend_param = list(at = c(-2, 2),
                                           tick_length = unit(0, "cm"),
                                           labels = c("Low", "High")),

               # Annotations
               top_annotation = ha,
               column_names_side = "top",
               column_names_rot = 45)

  # Draw the heatmap
  leg = Legend(col_fun = col_fun, title = "Average\nExpression", legend_height = unit(4, "cm"), grid_width = unit(1, "cm"))

  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
}
