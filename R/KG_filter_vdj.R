#' Filter V(D)J-seq data.
#'
#' This function streamlines the filtering process of V(D)J-seq data for in-depth BCR analysis.
#' @param vdj_data A data frame of V(D)J-seq data, or a list of data frames containing multiple V(D)J-seq.
#' @param multiple_samples Logical; whether or not vdj_data consists of a single data frame or a list of data frames.
#' @param remove_duplicate_IGH Logical; whhether heavy chain duplicates should be removed.
#' @param remove_duplicate_IGK Logical; whether kappa light chain duplicates should be removed.
#' @param remove_duplicate_IGL Logical; whether lambda light chain duplicates should be removed.
#' @param remove_no_c Logical; whether sequences with no annotated C gene should be removed.
#' @param remove_no_j Logical; whether sequences with no annotated J gene should be removed.
#' @param remove_no_v Logical; whether sequences with no annotated V gene should be removed
#' @param productive Logical; whether only productively recombined sequuences should be retained.
#' @param only_paired Logical; whether to only keep sequences that have a matching pair of heavy and light chains.
#' @param save_filtered Logical; whether sequences that are removed should be saved to an xlsx file.
#' @param print_numbers Logical; whether the number of sequences removed in each step should be printed in the console.
#' @export
#' @examples
#' KG_filter_vdj(vdj_dataset)
#'
KG_filter_vdj = function(vdj_data,
                         multiple_samples = FALSE,
                         remove_duplicate_IGH = TRUE,
                         remove_duplicate_IGK = TRUE,
                         remove_duplicate_IGL = TRUE,
                         remove_no_c = TRUE,
                         remove_no_j = TRUE,
                         remove_no_v = TRUE,
                         productive = TRUE,
                         only_paired = TRUE,
                         save_filtered = FALSE,
                         print_numbers = TRUE){

  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(openxlsx))

  if(multiple_samples){
    for(i in 1:length(vdj_data)){

      # Remove unproductive sequences
      if(productive){
        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], productive == "false"),
                               paste0(names(vdj_data)[i], "_unproductive_sequences.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of unproductive sequences in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], productive == "false")), "\n"))
        }

        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], productive == "true")
      }

      # Remove duplicate heavy chains
      if(remove_duplicate_IGH){
        multi_heavy = table(dplyr::filter(vdj_data[[i]], locus == "IGH")$cell_id)
        multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], cell_id %in% multi_heavy_cells),
                               paste0(names(vdj_data)[i], "_multiple_IGH_sequences.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of sequences with >1 IGH in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], cell_id %in% multi_heavy_cells)), "\n"))
        }

        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% multi_heavy_cells)
      }

      # Remove duplicate kappa light chains
      if(remove_duplicate_IGK){
        multi_igk = table(dplyr::filter(vdj_data[[i]], locus == "IGK")$cell_id)
        multi_igk_cells = names(multi_igk)[multi_igk > 1]

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], cell_id %in% multi_igk_cells),
                               paste0(names(vdj_data)[i], "_multiple_IGK_sequences.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of sequences with >1 IGK in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], cell_id %in% multi_igk_cells)), "\n"))
        }
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% multi_igk_cells)
      }

      # Remove duplicate lambda light chains
      if(remove_duplicate_IGL){
        multi_igl = table(dplyr::filter(vdj_data[[i]], locus == "IGL")$cell_id)
        multi_igl_cells = names(multi_igl)[multi_igl > 1]

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], cell_id %in% multi_igl_cells),
                               paste0(names(vdj_data)[i], "_multiple_IGL_sequences.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of sequences with >1 IGL in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], cell_id %in% multi_igl_cells)), "\n"))
        }
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% multi_igl_cells)
      }

      # Remove sequences with no annotated c_call
      if(remove_no_c){
        no_c = dplyr::filter(vdj_data[[i]], c_call == "")$cell_id

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], cell_id %in% no_c),
                               paste0(names(vdj_data)[i], "_no_c_call.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of sequences without c_call in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], cell_id %in% no_c)), "\n"))
        }
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% no_c)
      }

      # Remove sequences with no annotated v_call
      if(remove_no_v){
        no_v = dplyr::filter(vdj_data[[i]], v_call == "")$cell_id

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], cell_id %in% no_v),
                               paste0(names(vdj_data)[i], "_no_v_call.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of sequences without v_call in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], cell_id %in% no_v)), "\n"))
        }
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% no_v)
      }

      # Remove sequences with no annotated j_call
      if(remove_no_j){
        no_j = dplyr::filter(vdj_data[[i]], j_call == "")$cell_id

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], cell_id %in% no_j),
                               paste0(names(vdj_data)[i], "_no_j_call.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of sequences without j_call in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], cell_id %in% no_j)), "\n"))
        }
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% no_j)
      }

      # Remove single heavy/light chain sequences
      if(only_paired){
        heavy_cells = dplyr::filter(vdj_data[[i]], locus == "IGH")$cell_id
        light_cells = dplyr::filter(vdj_data[[i]], locus == "IGL" | locus == "IGK")$cell_id
        paired = intersect(heavy_cells, light_cells)

        if(save_filtered){
          openxlsx::write.xlsx(dplyr::filter(vdj_data[[i]], !cell_id %in% paired),
                               paste0(names(vdj_data)[i], "_unpaired_sequences.xlsx"))
        }
        if(print_numbers){
          cat(paste0("Number of unpaired sequences in ", names(vdj_data)[i], ": ", nrow(dplyr::filter(vdj_data[[i]], !cell_id %in% paired)), "\n"))
        }
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], cell_id %in% paired)
      }
    }
  } else{

    # Remove unproductive sequences
    if(productive){
      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, productive == "false"),
                             paste0(deparse(vdj_data), "_unproductive_sequences.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of unproductive sequences: ", nrow(dplyr::filter(vdj_data, productive == "false")), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, productive == "true")
    }
    # Remove duplicate heavy chains
    if(remove_duplicate_IGH){
      multi_heavy = table(dplyr::filter(vdj_data, locus == "IGH")$cell_id)
      multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, cell_id %in% multi_heavy_cells),
                             paste0(deparse(vdj_data), "_multiple_IGH_sequences.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of sequences with >1 IGH: ", nrow(dplyr::filter(vdj_data, cell_id %in% multi_heavy_cells)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% multi_heavy_cells)
    }

    # Remove duplicate kappa light chains
    if(remove_duplicate_IGK){
      multi_igk = table(dplyr::filter(vdj_data, locus == "IGK")$cell_id)
      multi_igk_cells = names(multi_igk)[multi_igk > 1]

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, cell_id %in% multi_igk_cells),
                             paste0(deparse(vdj_data), "_multiple_IGK_sequences.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of sequences with >1 IGK: ", nrow(dplyr::filter(vdj_data, cell_id %in% multi_igk_cells)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% multi_igk_cells)
    }

    # Remove duplicate lambda light chains
    if(remove_duplicate_IGL){
      multi_igl = table(dplyr::filter(vdj_data, locus == "IGL")$cell_id)
      multi_igl_cells = names(multi_igl)[multi_igl > 1]

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, cell_id %in% multi_igl_cells),
                             paste0(deparse(vdj_data), "_multiple_IGL_sequences.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of sequences with >1 IGL: ", nrow(dplyr::filter(vdj_data, cell_id %in% multi_igl_cells)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% multi_igl_cells)
    }

    # Remove sequences with no annotated c_call
    if(remove_no_c){
      no_c = dplyr::filter(vdj_data, c_call == "")$cell_id

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, cell_id %in% no_c),
                             paste0(deparse(vdj_data), "_no_c_call.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of sequences without c_call in ", deparse(vdj_data), ": ", nrow(dplyr::filter(vdj_data, cell_id %in% no_c)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% no_c)
    }

    # Remove sequences with no annotated v_call
    if(remove_no_v){
      no_v = dplyr::filter(vdj_data, v_call == "")$cell_id

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, cell_id %in% no_v),
                             paste0(deparse(vdj_data), "_no_v_call.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of sequences without v_call in ", deparse(vdj_data), ": ", nrow(dplyr::filter(vdj_data, cell_id %in% no_v)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% no_v)
    }

    # Remove sequences with no annotated j_call
    if(remove_no_j){
      no_j = dplyr::filter(vdj_data, j_call == "")$cell_id

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, cell_id %in% no_j),
                             paste0(deparse(vdj_data), "_no_j_call.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of sequences without j_call in ", deparse(vdj_data), ": ", nrow(dplyr::filter(vdj_data, cell_id %in% no_j)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% no_j)
    }

    # Remove single heavy/light chain sequences
    if(only_paired){
      heavy_cells = dplyr::filter(vdj_data, locus == "IGH")$cell_id
      light_cells = dplyr::filter(vdj_data, locus == "IGL" | locus == "IGK")$cell_id
      paired = intersect(heavy_cells, light_cells)

      if(save_filtered){
        openxlsx::write.xlsx(dplyr::filter(vdj_data, !cell_id %in% paired),
                             paste0(deparse(vdj_data), "_unpaired_sequences.xlsx"))
      }
      if(print_numbers){
        cat(paste0("Number of unpaired sequences: ", nrow(dplyr::filter(vdj_data, !cell_id %in% paired)), "\n"))
      }
      vdj_data = dplyr::filter(vdj_data, cell_id %in% paired)
    }
    gc()
  }
  return(vdj_data)
}
