#' Calling clones in BCR data
#'
#' This function is a wrapper for calling clonotypes in BCR data using the dowser and scoper packages.
#' @param bcr_data A data frame containing the BCR data in AIRR format.
#' @param threshold Threshold used for hierarchical clustering of sequences. Default is 0.15, i.e., 85% similarity threhsold.
#' @param reference_dir Path to directory with IMGT VDJ references to use.
#' @export
#' @examples
#' KG_clones(bcr_dataset)
KG_clones = function(bcr_data,
                     threshold = 0.15,
                     reference_dir = NULL){
  suppressPackageStartupMessages(c(require(scoper),
                                   require(dowser),
                                   require(shazam),
                                   require(dplyr)))

  # Cluster BCR sequences
  cat("Clustering BCR sequences")
  res = scoper::hierarchicalClones(bcr_data,
                                   cell_id = "cell_id",
                                   threshold = threshold,
                                   only_heavy = FALSE,
                                   split_light = TRUE,
                                   summarize_clones = FALSE)

  # Read in references
  references = dowser::readIMGT(dir = reference_dir)

  igh = dowser::createGermlines(filter(res, locus == "IGH"), references)
  igk = dowser::createGermlines(filter(res, locus == "IGK"), references)
  igl = dowser::createGermlines(filter(res, locus == "IGL"), references)

  # Calculate SHM frequency
  ## Heavy chain
  cat("Calculating SHM frequency in...")
  cat("heavy chains")
  igh_data = shazam::observedMutations(igh,
                                       sequenceColumn = "sequence_alignment",
                                       germlineColumn = "germline_alignment_d_mask",
                                       regionDefinition = IMGT_V,
                                       frequency = TRUE,
                                       combine = TRUE,
                                       nproc = 1)
    ## Kappa chain
  cat("kappa chains")
  igk_data = shazam::observedMutations(igk,
                                       sequenceColumn = "sequence_alignment",
                                       germlineColumn = "germline_alignment_d_mask",
                                       regionDefinition = IMGT_V,
                                       frequency = TRUE,
                                       combine = TRUE,
                                       nproc = 1)
  ## Lambda chain
  cat("lambda chains")
  igl_data = shazam::observedMutations(igl,
                                       sequenceColumn = "sequence_alignment",
                                       germlineColumn = "germline_alignment_d_mask",
                                       regionDefinition = IMGT_V,
                                       frequency = TRUE,
                                       combine = TRUE,
                                       nproc = 1)

  # Combine data
  clone_data = rbind(igh_data, igl_data, igk_data)
cat("Done")
  return(clone_data)
}
