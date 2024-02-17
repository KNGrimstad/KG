#' Add suffix to list of barcodes
#' This function adds a suffix to all cells in the specified input list.
#' @param cells A list of cell names.
#' @param suffix The suffix to be added to the cell names
#' @export
#' @examples
#' KG_add_suffix(list_of_cells, "_4")
KG_add_suffix = function(cells, suffix){
  paste0(cells, suffix)
}
