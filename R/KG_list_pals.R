#' List all available color palettes in KG.
#'
#' This function lists the names of all color palettes available in this package.
#'
#' @return A character vector of palette names.
#' @export
#' @examples
#' KG_list_palettes()
KG_list_pals = function() {

  pal_functions <- ls(envir = asNamespace("KG"), pattern = "^KG_.*pal$")

  gsub("^KG_", "", pal_functions)
}
