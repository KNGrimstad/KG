#' A blinding palette!
#'
#' This function gnereates an illuminating color palette ranging from black, to gold, to white.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_headlight_pal(50)
KG_headlight_pal = function(n){
  suppressPackageStartupMessages(require(grDevices))
  grDevices::colorRampPalette(c("black", "gold", "white"))(n)
}
