#' A prideful rainbow palette!
#'
#' This function generates an inclusive palette for everyone.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_flame_pal(50)
KG_rainbow_pal = function(n){
  suppressPackageStartupMessages(require(grDevices))
  grDevices::colorRampPalette(c("#FF0000", "#FFA510", "#FFFF04", "#00FF04", "#1E90F2", "#4B369D", "#70369E"))(n)
}
