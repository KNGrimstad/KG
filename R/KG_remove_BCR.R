#' Remove BCR genes
#'
#' This function removes all B-cell receptor genes from a Seurat object.
#' @param seurat_object A Seurat object
#' @export
#' @examples
#' KG_remove_BCR(B_cell_dataset)
KG_remove_BCR = function(seurat_object){
  seurat_object = seurat_object[!grepl("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]|^AC233755.1$",
                                       rownames(seurat_object)),]
  return(seurat_object)
}
