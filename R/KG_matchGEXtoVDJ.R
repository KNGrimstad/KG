#' Match V(D)J-seq sequences to GEX data
#'
#' This function filters V(D)J-seq data based on sequences present in paired GEX data. Sequences that are not present in the GEX data are removed.
#' @param vdj_data An AIRR data object, or a list of AIRR data objects.
#' @param gex_data A Seurat object, or a list of Seurat objects
#' @param multiple_samples Logical; whether or not vdj_data consists of a single data frame or a list of data frames. NOTE: if set to TRUE, the function expects gex_data is also a list of matching Seurat objects. Matching objects should appear in the same order in both GEX and V(D)J-seq lists.
#' @param save_filtered Logical; whether sequences that are removed should be saved to an xlsx file.
#' @param print_numbers Logical; whether the number of sequences removed in each step should be printed in the console.
#' @export
#' @examples
#' KG_matchGEXtoVDJ(seurat_object, vdj_dataset)
#'
KG_matchGEXtoVDJ = function(gex_data,
                            vdj_data,
                            multiple_samples = FALSE,
                            save_filtered = FALSE,
                            print_numbers = TRUE){
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(openxlsx))
  suppressPackageStartupMessages(require(Seurat))

  if(multiple_samples){
    if(is.null(names(vdj_data))){
      names(vdj_data) = paste0("vdj_data", seq(length(vdj_data)))
    }
    # If vdj_data is a list of data frames
    for(i in 1:length(vdj_data)){
      nseq_vdj = length(unique(vdj_data[[i]][["cell_id"]]))
      nseq_gex = ncol(gex_data[[i]])
      matching_seqs = intersect(Seurat::Cells(gex_data[[i]]), vdj_data[[i]][["cell_id"]])
      percent_match = round(length(matching_seqs) / nseq_gex * 100, 1)
      if(print_numbers){
        cat(paste0("Total number of V(D)J sequences in ", names(vdj_data)[i], ": ", nseq_vdj, "\n"))
        cat(paste0("Number of matching sequences in ", names(vdj_data)[i], ": ", length(matching_seqs), " (", percent_match, "% of GEX matching V(D)J)\n"))
      }
      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], !cell_id %in% matching_seqs),
                             paste0(names(vdj_data)[i], "_sequences_not_matching_GEX.xlsx"))
      }

      # Remove sequences not present in GEX data
      vdj_data[[i]] = dplyr::filter(vdj_data[[i]], cell_id %in% matching_seqs)
    }
  } else{
    # If vdj_data is a single data frame
    nseq_vdj = length(unique(vdj_data[["cell_id"]]))
    nseq_gex = ncol(gex_data)
    matching_seqs = intersect(Seurat::Cells(gex_data), vdj_data[["cell_id"]])
    percent_match = round(length(matching_seqs) / nseq_gex * 100, 1)
    if(print_numbers){
      cat(paste0("Total number of V(D)J sequences: ", length(nseq_vdj), "\n"))
      cat(paste0("Number of matching sequences: ", length(matching_seqs), " (", percent_match, "% of GEX matching V(D)J)\n"))
    }
    if(save_filtered){
      openxlsx::write.xlsx(dplyr::filter(vdj_data, !cell_id %in% matching_seqs),
                           paste0(substitute(vdj_data), "_sequences_not_matching_GEX.xlsx"))
    }
  }
  gc()
  return(vdj_data)
}
