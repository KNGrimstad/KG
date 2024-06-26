#' Hierarchical clustering
#'
#' This function is a wrapper for hierarchical clustering of Seurat clusters.
#' @param seurat_object A Seurat object
#' @param genes A list of genes to base the clustering on. Default is the variable features.
#' @param assay The name of the Seurat object assay to use.
#' @param title Plot title.
#' @param ident_col Identity column to use.Default is 'seurat_clusters'.
#' @param xlim Determines x-axis limits.
#' @export
#' @examples
#' KG_hc(B_cell_dataset)
KG_hc = function(seurat_object,
                 genes = NULL,
                 assay = "RNA",
                 title = NULL,
                 ident_col = "seurat_clusters",
                 xlim = c(NA, 200)){
  suppressPackageStartupMessages({
    require(Seurat)
    require(SeuratObject)
    require(ggtree)
    require(ggplot2)
  })

  temp = seurat_object
  Seurat::Idents(temp) = temp[[ident_col, drop = TRUE]]

  # Construct a cluster tree
  if(is.null(genes)){
    genes = Seurat::VariableFeatures(seurat_object)
  }
  cat("Building cluster tree")
  tree = Seurat::BuildClusterTree(seurat_object,
                                  assay = assay,
                                  features = genes,
                                  verbose = FALSE)

  phytree = SeuratObject::Tool(object = tree,
                 slot = "BuildClusterTree")

  if(is.null(title)){
    dendro = ggtree::ggtree(phytree) +
      ggtree::geom_tiplab() +
      ggtree::theme_tree() +
      ggtree::xlim(xlim)
  } else {
    dendro = ggtree::ggtree(phytree) +
      ggtree::geom_tiplab() +
      ggtree::theme_tree() +
      ggtree::xlim(xlim) +
      ggplot2::ggtitle(title)
  }
  return(dendro)
}
