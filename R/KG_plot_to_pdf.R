#' Print plots trees to PDF.
#'
#' This function prints a list of plots to a PDF file.
#' @param plots A list of plots.
#' @param ncol The number of columns per page.
#' @param nrow The number of rows per page.
#' @param heading The heading to be printed on each page.
#' @param padding Distance between individual plots.
#' @param file The name of the PDF file.
#' @param pdf_size Size of the PDF.
#' @param portrait Wheter or not to use portrait layout (if FALSE, will use horizontal layout).
#' @export
#' @examples
#' KG_trees_to_pdf(plots, file = "plot_file")
KG_plot_to_pdf = function(plots,
                           ncol = 3,
                           nrow = 5,
                           heading = NULL,
                           padding = NULL,
                           file = NULL,
                           pdf_size = a4,
                           portrait = TRUE){

  suppressPackageStartupMessages({
    require(stats)
    require(grDevices)
    require(grid)})

  # Formats
  a4 = c(8.27, 11.69)
  a5 = c(5.83, 8.27)
  a6 = c(4.13, 5.83)

  # Initiate PDF
  if(portrait){
    grDevices::pdf(paste(file, ".pdf", sep = ""),
                   width = pdf_size[1], height = pdf_size[2])
  } else{
    grDevices::pdf(paste(file, ".pdf", sep = ""),
                   width = pdf_size[2], height = pdf_size[1])
  }

  # Print plots
  gridExtra::marrangeGrob(grobs = plots,
                          ncol = ncol,
                          nrow = nrow,
                          top = heading,
                          padding = if(is.null(padding)) grid::unit(1, "lines") else padding,
                          pages = length(plots) %/% (ncol * nrow) + ifelse(length(plots) %% (ncol * nrow) > 0, 1, 0)) %>%
    grid::grid.draw()
  grDevices::dev.off()
}
