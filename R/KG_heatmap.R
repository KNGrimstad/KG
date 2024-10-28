#' Heatmap of relative gene expression
#'
#' This functions draws a heatmap displaying the relative cluster-average expression of specific genes.
#' @param seurat_object A Seurat object.
#' @param genes A vector of genes for which to plot the expression.
#' @param assay The name of the assay to extract expression values from.
#' @param group_by Group variable, if NULL, uses the active identity.
#' @param show_genes Whether or not to show the gene names on the y-axis.
#' @param plot_only_once Logical; whether a gene should be plotted only once or each time it appears in the genes vector.
#' @param expression_cols A vector of colors used for plotting the relative expression.
#' @param cols A vector of colors to use for the group bars. If NULL, default colors are used.
#' @param heading Whether or not to include a heaing to the heatmap.
#' @param title The title, if heading = TRUE.
#' @param column_title_size Font size of title text.
#' @param group_label_size Font size of group labels.
#' @param row_label_size Font size of rownames (genes). If NA, will set dynamic font size.
#' @param group_border_col The color of the border for the group label boxes.
#' @param group_border_lwd Line width of the border for the group label boxes.
#' @param border Logical; whether borders should be drawn.
#' @param border_col_major The color for the outer borders of the heatmap.
#' @param border_lwd_major Line width for the outer borders of the heatmap.
#' @param border_col_minor The color for the inner borders of the heatmap (the tiles).
#' @param border_lwd_minor Line width for the inner borders of the heatmap (the tiles).
#' @param cluster_rows Whether or not genes should be clustered.
#' @param cluster_columns Whether or not groups should be clustered.
#' @param row_splits A value, or vector of values, detailing the group number and the number of genes in the group, to split rows between groups, e.g., c(rep(1, 3), rep(2, 4)) adds gaps between the third and the fourth gene, assuming there are seven genes in total.
#' @param row_gap_size The size o the row gap in mm.
#' @param column_splits Same as split_row_at, but for columns.
#' @param raster Whether or not the heatmap should be rasterized
#' @export
#' @examples
#' KG_heatmap(B_cell_dataset, genes = c("CD19", "CD27", "HLA-DRB1"))
KG_heatmap = function (seurat_object,
                       genes,
                       assay = "RNA",
                       group_by = NULL,
                       show_genes = TRUE,
                       plot_only_once = TRUE,
                       expression_cols = c("#FF00FF", "black", "#FFFF00"),
                       cols = NULL,
                       heading = TRUE,
                       title = NA,
                       column_title_size = 15,
                       group_label_size = 15,
                       row_label_size = NULL,
                       group_border_col = "black",
                       group_border_lwd = 1,
                       border = TRUE,
                       border_col_major = "black",
                       border_lwd_major = 1,
                       border_col_minor = "black",
                       border_lwd_minor = 1,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       row_splits = NULL,
                       row_gap_size = 2,
                       column_splits = NULL,
                       raster = FALSE){

  suppressPackageStartupMessages({
    require(ComplexHeatmap)
    require(scales)
    require(circlize)
    require(Seurat)
    require(grid)
  })
  col_fun = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = expression_cols)
  if (!is.null(group_by)) {
    Seurat::Idents(seurat_object) = factor(seurat_object[[group_by,
                                                          drop = TRUE]])
  }
  if(plot_only_once){
    if(is.null(genes)){
      stop("No genes included.")
    }
    avg_matrix = t(scale(t(as.matrix(Seurat::AverageExpression(seurat_object,
                                                               features = genes)[[assay]]))))
    annotations_df = data.frame(factor(levels(Seurat::Idents(seurat_object))))
    names(annotations_df) = "Cluster"
    rownames(annotations_df) = levels(Seurat::Idents(seurat_object))

    if (is.null(cols)) {
      cols = (scales::hue_pal())(length(unique(annotations_df$Cluster)))
    }
    names(cols) = unique(annotations_df$Cluster)

    ha = ComplexHeatmap::HeatmapAnnotation(df = data.frame(Cluster = annotations_df$Cluster),
                                           col = list(Cluster = cols),
                                           which = "col",
                                           show_annotation_name = FALSE,
                                           border = TRUE,
                                           gp = grid::gpar(col = group_border_col,
                                                           lwd = group_border_lwd),
                                           show_legend = FALSE)
    if(is.null(row_label_size)){
      dynamic_font = max(5, min(15, 20 - length(genes)/5))
    }

    ht = ComplexHeatmap::Heatmap(avg_matrix,
                                 border = border,
                                 border_gp = grid::gpar(col = border_col_major,
                                                        lty = 1, lwd = border_lwd_major),
                                 rect_gp = grid::gpar(col = border_col_minor,
                                                      lwd = border_lwd_minor),
                                 col = col_fun,
                                 use_raster = raster,

                                 # Rows
                                 cluster_rows = cluster_rows,
                                 show_row_names = show_genes,
                                 row_title = NULL,
                                 row_dend_reorder = TRUE,
                                 row_split = row_splits,
                                 row_names_side = "left",
                                 row_names_gp = grid::gpar(fontsize = if (is.null(row_label_size)) dynamic_font else row_label_size),
                                 row_gap = unit(row_gap_size, "mm"),
                                 #gap = unit(split_row_at),

                                 # Columns
                                 cluster_columns = cluster_columns,
                                 column_title_gp = grid::gpar(fontsize = column_title_size),
                                 column_names_gp = grid::gpar(fontsize = group_label_size),
                                 show_column_names = TRUE,
                                 column_title = if (!is.na(title) && heading) title else NULL,
                                 column_dend_reorder = TRUE,
                                 column_split = column_splits,

                                 # Legend
                                 name = "Average\nExpression",
                                 heatmap_legend_param = list(at = c(-2, 2),
                                                             tick_length = grid::unit(0, "cm"),
                                                             labels = c("Low", "High")),
                                 top_annotation = ha, column_names_side = "top",
                                 column_names_rot = 45)

    leg = ComplexHeatmap::Legend(col_fun = col_fun,
                                 title = "Average\nExpression",
                                 legend_height = grid::unit(4, "cm"),
                                 grid_width = grid::unit(1, "cm"))

    ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                         annotation_legend_side = "right")

  } else{

    if(is.null(genes)){
      stop("No genes included")
    }

    # Create unique identifiers for genes
    gene_counts = table(genes)
    indexed_genes = genes
    for(gene in names(gene_counts)){
      if(gene_counts[gene] > 1){
        indices = which(genes == gene)
        indexed_genes[indices] = paste0(gene, "_", seq_along(indices))
      }
    }
    avg_matrix = t(scale(t(as.matrix(Seurat::AverageExpression(seurat_object,
                                                               features = genes)[[assay]]))))

    new_avg_matrix = matrix(nrow = length(genes), ncol = ncol(avg_matrix)) # Create a new matrix
    rownames(new_avg_matrix) = indexed_genes
    colnames(new_avg_matrix) = colnames(avg_matrix)

    for(i in seq_along(indexed_genes)){
      base_gene = gsub("_\\d+", "", indexed_genes[i])
      new_avg_matrix[i, ] = avg_matrix[base_gene, ]
    }

    avg_matrix = new_avg_matrix

    annotations_df = data.frame(factor(levels(Seurat::Idents(seurat_object))))
    names(annotations_df) = "Cluster"
    rownames(annotations_df) = levels(Seurat::Idents(seurat_object))

    if (is.null(cols)) {
      cols = (scales::hue_pal())(length(unique(annotations_df$Cluster)))
    }
    names(cols) = unique(annotations_df$Cluster)

    ha = ComplexHeatmap::HeatmapAnnotation(df = data.frame(Cluster = annotations_df$Cluster),
                                           col = list(Cluster = cols),
                                           which = "col",
                                           show_annotation_name = FALSE,
                                           border = TRUE,
                                           gp = grid::gpar(col = group_border_col,
                                                           lwd = group_border_lwd),
                                           show_legend = FALSE)

    if(is.null(row_label_size)){
      dynamic_font = max(5, min(15, 20 - length(genes)/5))
    }

    ht = ComplexHeatmap::Heatmap(avg_matrix,
                                 border = border,
                                 border_gp = grid::gpar(col = border_col_major,
                                                        lty = 1, lwd = border_lwd_major),
                                 rect_gp = grid::gpar(col = border_col_minor, lwd = border_lwd_minor),
                                 col = col_fun,
                                 use_raster = raster,

                                 # Rows
                                 cluster_rows = cluster_rows,
                                 show_row_names = show_genes,
                                 row_title = NULL,
                                 row_dend_reorder = TRUE,
                                 row_split = row_splits,
                                 row_names_side = "left",
                                 row_names_gp = grid::gpar(fontsize = if (is.null(row_label_size)) dynamic_font else row_label_size),
                                 row_gap = unit(row_gap_size, "mm"),
                                 row_labels = genes,
                                 #gap = unit(split_row_at),

                                 # Columns
                                 cluster_columns = cluster_columns,
                                 column_title_gp = grid::gpar(fontsize = column_title_size),
                                 column_names_gp = grid::gpar(fontsize = group_label_size),
                                 show_column_names = TRUE,
                                 column_title = if (!is.na(title) && heading) title else NULL,
                                 column_dend_reorder = TRUE,
                                 column_split = column_splits,

                                 # Legend
                                 name = "Average\nExpression",
                                 heatmap_legend_param = list(at = c(-2, 2),
                                                             tick_length = grid::unit(0, "cm"),
                                                             labels = c("Low", "High")),
                                 top_annotation = ha, column_names_side = "top",
                                 column_names_rot = 45)

    leg = ComplexHeatmap::Legend(col_fun = col_fun,
                                 title = "Average\nExpression",
                                 legend_height = grid::unit(4, "cm"),
                                 grid_width = grid::unit(1, "cm"))

    ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                         annotation_legend_side = "right")
  }
}
