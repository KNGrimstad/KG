#' A flaming hot color palette!
#'
#' This function generates a fiery color palette, ranging from black, through red and orange, to gold.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_flame_pal(50)
KG_flame_pal = function(n){
  suppressPackageStartupMessages(require(grDevices))
  grDevices::colorRampPalette(c("black", "red3", "orange", "gold"))(n)
}
