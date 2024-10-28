#' A slimy color palette for ninja mutatns!
#'
#' This function generates a radioactive color palette, ranging from black through successively brighter shades of green.
#' @param n The number of colors to generate.
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' KG_mutant_pal(50)
KG_mutant_pal = function(n){
  mutant = colorRampPalette(c("black", "darkgreen", "limegreen", "springgreen"))(n)
}
