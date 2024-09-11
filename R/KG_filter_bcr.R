#' Filter AIRR object BCR data
#'
#' This function filters BCR data output from IgBLAST.
#' @param bcr_data A data frame containing the BCR data in AIRR format.
#' @param productive Whether unproductive transcripts should be removed.
#' @param single_IgH Whether cells with multiple heavy chains should be removed.
#' @param c_annotated Whether cells with no annotated C gene should be removed.
#' @param hl_pair Whether cells with no paired heavy/light chain should be removed.
#' @export
#' @examples
#' KG_filter_bcr(BCR_object)
KG_filter_bcr = function(bcr_data,
                         productive = TRUE,
                         single_IgH = TRUE,
                         c_annotated = TRUE,
                         hl_pair = TRUE){

  suppressPackageStartupMessages({
    require(dplyr)
    require(stats)})

  warning("KG_filter_bcr is deprecated, please use KG_filter_vdj instead.")

  # Count number of sequences before filtering
  base_cells = nrow(bcr_data)

  # Remove unproductive chains
  if(productive){
    message("Removing unproductive transcripts")
    bcr_data = bcr_data %>%
      dplyr::filter(productive)
    message(paste(base_cells - nrow(bcr_data), " sequences were removed.", sep = ""))
  }

  # Remove cells with multiple heavy chains
  if(single_IgH){
    message("Removing cells with multiple heavy chains.")
    base_cells1 = nrow(bcr_data)
    multi_h = table(dplyr::filter(bcr_data, locus == "IGH")$cell_id)
    multi_h_cells = names(multi_h)[multi_h >1]
    bcr_data = dplyr::filter(bcr_data, !cell_id %in% multi_h_cells)
    message(paste(base_cells1 - nrow(bcr_data), " sequences were removed.", sep = ""))
  }

  # Remove cells with no annotated C gene
  if(c_annotated){
    base_cells2 = nrow(bcr_data)
    message("Removing cells with no annotated C gene.")
    bcr_data = bcr_data[stats::complete.cases(bcr_data$c_call),]
    message(paste(base_cells2 - nrow(bcr_data), " sequences were removed.", sep = ""))
  }

  # Remove unpaired sequences
  if(hl_pair){
    message("Removing cells with unpaired IgH/IgL chains.")
    base_cells3 = nrow(bcr_data)
    h_cells = dplyr::filter(bcr_data, locus == "IGH")$cell_id
    l_cells = dplyr::filter(bcr_data, locus == "IGK" | locus == "IGL")$cell_id
    no_heavy = l_cells[which(!l_cells %in% h_cells)]
    no_light = h_cells[which(!h_cells %in% l_cells)]
    bcr_data = dplyr::filter(bcr_data, !cell_id %in% c(no_heavy, no_light))
    message(paste(base_cells3 - nrow(bcr_data), " sequences were removed.", sep = ""))
  }
  message(paste("Filtering complete. ", ba