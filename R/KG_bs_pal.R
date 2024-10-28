#' A color palette that swings both ways.
#'
#' This function generates a color palette ranging from pink to purple to dark blue.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_bs_pal(50)
KG_bs_pal = function(n){
  suppressPackageStartupMessages(require(grDevices))
  grDevices::colorRampPalette(c(hcl(330, 75, 50),
                     hcl(2450, 75, 35),
                     hcl(240, 75, 15)),
                   bias = 1)(n)
}
