#' Clonotype analysis
#'
#' This function streamlines the functions from the Immcantation framework to define hierarchical clones.
#' @param vdj_data An AIRR data object.
#' @param cell_id Name of the column containing the cell identifiers/barcodes.
#' @param only_heavy Logical; whether to only use IGH sequences for grouping.
#' @param split_light Logical; whether to split clones for light chains.
#' @param summarize_clones Logical; same as in scoper::hierarchicalClones.
#' @param fields Vector of additional grouping columns. Sequences with different values will be considered as separate clones. Useful for separation by sample in integrated datasets.
#' @export
#' @examples
#' KG_prep_lones(vdj_data)
#'
KG_prep_clones = function(vdj_data,
                          cell_id = "cell_id",
                          only_heavy = TRUE,
                          split_light = TRUE,
                          summarize_clones = FALSE,
                          fields = NULL){

  suppressPackageStartupMessages({
    require(shazam)
    require(dowser)
    require(scoper)
    require(dplyr)
    require(ggplot2)})

  # Distance to nearest
  dist_nearest = shazam::distToNearest(dplyr::filter(vdj_data, locus == "IGH"))

  ggplot2::ggplot(dplyr::filter(dist_nearest, !is.na(dist_nearest)),
                  ggplot2::aes(x = dist_nearest)) +
    ggplot2::geom_histogram(color = "white", binwidth = 0.02) +
    ggplot2::labs(x = "Hamming distance", y = "Count") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 16))

  thresh = readline(prompt = "Select distance threshold for clonal grouping: \n")

  results = scoper::hierarchicalClones(vdj_data,
                                       cell_id = cell_id,
                                       threshold = thresh,
                                       only_heavy = only_heavy,
                                       split_light = split_light,
                                       summarize_clones = summarize_clones,
                                       fields = fields,
                                       verbose = TRUE)

  return(results)
}
