#' Calling clones in BCR data
#' This function calls clonotypes in BCR data using the dowser and scoper packages.
#' @param bcr_data A data frame containing the BCR data in AIRR format.
#' @param threshold Threshold used for hierarchical clustering of sequences. Default is 0.15, i.e., 85% similarity threhsold.
#' @param reference_dir Path to directory with IMGT VDJ references to use.
#' @export
#' @examples
#' KG_clones(bcr_dataset)
KG_clones = function(bcr_data,
                     threshold = 0.15,
                     reference_dir = NULL){
  require(scoper)
  require(dowser)
  # Cluster BCR sequences
  message("Clustering BCR sequences")
  res = hierarchicalClones(bcr_data,
                           cell_id = "cell_id",
                           threshold = threshold,
                           only_heavy = FALSE,
                           split_light = TRUE,
                           summarize_clones = FALSE)

  # Read in references
  references = readIMGT(dir = reference_dir)

  igh = createGermlines(filter(res, locus == "IGH"), references)
  igk = createGermlines(filter(res, locus == "IGK"), references)
  igl = createGermlines(filter(res, locus == "IGL"), references)

  # Calculate SHM frequency
  ## Heavy chain
  message("Calculating SHM frequency in\nheavy chains")
  igh_data = observedMutations(igh,
                               sequenceColumn = "sequence_alignment",
                               germlineColumn = "germline_alignment_d_mask",
                               regionDefinition = IMGT_V,
                               frequency = TRUE,
                               combine = TRUE,
                               nproc = 1)
    ## Kappa chain
  message("kappa chains")
  igk_data = observedMutations(igk,
                             sequenceColumn = "sequence_alignment",
                             germlineColumn = "germline_alignment_d_mask",
                             regionDefinition = IMGT_V,
                             frequency = TRUE,
                             combine = TRUE,
                             nproc = 1)
  ## Lambda chain
  message("lambda chains")
  igl_data = observedMutations(igl,
                             sequenceColumn = "sequence_alignment",
                             germlineColumn = "germline_alignment_d_mask",
                             regionDefinition = IMGT_V,
                             frequency = TRUE,
                             combine = TRUE,
                             nproc = 1)

  # Combine data
  clone_data = rbind(igh_data, igl_data, igk_data)

  return(clone_data)
}
