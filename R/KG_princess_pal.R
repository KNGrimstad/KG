#' A sickeningly sweet bubblegum palette!
#'
#' This function generates the cotton candy palette of every princess' dream.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_flame_pal(50)
KG_princess_pal = function(n){
  suppressPackageStartupMessages(require(grDevices))
  grDevices::colorRampPalette(c("#FFFFFF", "#FFE1FF", "#FFC0CB", "#FF82AB", "#FF6EB4" , "#FF1493"))(n)
}
