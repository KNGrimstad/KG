#' Remove suffix from list of barcodes
#'
#' This function removes a suffix from all cells in the specified input list.
#' @param cells A list of cell names from a Seurat object.
#' @param suffix The suffix to be removed from cell names.
#' @export
#' @examples
#' KG_remove_suffix(list_of_cells, "_4")
KG_remove_suffix = function(cells, suffix){
  result = sub(paste0(suffix, "$"), "", cells)
  return(result)
}
