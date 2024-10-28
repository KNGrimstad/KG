#' An eery color palette!
#'
#' This function generates a color palette from black to successively brighter shades of blue.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_ghost_pal(50)
KG_ghost_pal = function(n){
  suppressPackageStartupMessages(require(grDevices))
  grDevices::colorRampPalette(c("black", "darkblue", "dodgerblue", "lightblue"))(n)
}
