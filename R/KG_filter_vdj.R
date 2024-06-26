#' Filter V(D)J-seq data.
#'
#' This function streamlines the filtering process of V(D)J-seq data for in-depth BCR analysis.
#' @param vdj_data An AIRR data object, or a list of AIRR data objects.
#' @param multiple_samples Logical; whether or not vdj_data consists of a single or multiple objects.
#' @param remove_duplicate_IGH Logical; whhether heavy chain duplicates should be removed.
#' @param remove_duplicate_IGK Logical; whether kappa light chain duplicates should be removed.
#' @param remove_duplicate_IGL Logical; whether lambda light chain duplicates should be removed.
#' @param remove_no_c Logical; whether sequences with no annotated C gene should be removed.
#' @param remove_no_j Logical; whether sequences with no annotated J gene should be removed.
#' @param remove_no_v Logical; whether sequences with no annotated V gene should be removed
#' @param productive Logical; whether only productively recombined sequuences should be retained.
#' @param only_paired Logical; whether to only keep sequences that have a matching pair of heavy and light chains.
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
                         only_paired = TRUE){

  suppressPackageStartupMessages(require(dplyr))

  if(multiple_samples){
    for(i in 1:length(vdj_data)){

      # Remove unproductive sequences
      if(productive){
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], productive == "true")
      }

      # Remove duplicate heavy chains
      if(remove_duplicate_IGH){
        mulit_heavy = table(dplyr::filter(vdj_data[[i]], locus == "IGH")$cell_id)
        multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% multi_heavy_cells)
      }

      # Remove duplicate kappa light chains
      if(remove_duplicate_IGK){
        multi_igk = table(dplyr::filter(vdj_data[[i]], locus == "IGK")$cell_id)
        multi_igk_cells = names(multi_igk)[multi_igk > 1]
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% multi_igk_cells)
      }

      # Remove duplicate lambda light chains
      if(remove_duplicate_IGL){
        multi_igl = table(dplyr::filter(vdj_data[[i]], locus == "IGL")$cell_id)
        multi_igl_cells = names(multi_igl)[multi_igl > 1]
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], !cell_id %in% multi_igl_cells)
      }

      # Remove single heavy/light chain sequences
      if(only_paired){
        heavy_cells = dplyr::filter(vdj_data[[i]], locus == "IGH")$cell_id
        light_cells = dplyr::filter(vdj_data[[i]], locus == "IGL" | locus == "IGK")$cell_id
        paired = intersect(heavy_cells, light_cells)
        vdj_data[[i]] = dplyr::filter(vdj_data[[i]], cell_id %in% paired)
      }
    }
  } else{

    # Remove unproductive sequences
    if(productive){
      vdj_data = dplyr::filter(vdj_data, productive == "true")
    }

    # Remove duplicate heavy chains
    if(remove_duplicate_IGH){
      mulit_heavy = table(dplyr::filter(vdj_data, locus == "IGH")$cell_id)
      multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% multi_heavy_cells)
    }

    # Remove duplicate kappa light chains
    if(remove_duplicate_IGK){
      multi_igk = table(dplyr::filter(vdj_data, locus == "IGK")$cell_id)
      multi_igk_cells = names(multi_igk)[multi_igk > 1]
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% multi_igk_cells)
    }

    # Remove duplicate lambda light chains
    if(remove_duplicate_IGL){
      multi_igl = table(dplyr::filter(vdj_data, locus == "IGL")$cell_id)
      multi_igl_cells = names(multi_igl)[multi_igl > 1]
      vdj_data = dplyr::filter(vdj_data, !cell_id %in% multi_igl_cells)
    }

    # Remove single heavy/light chain sequences
    if(only_paired){
      heavy_cells = dplyr::filter(vdj_data, locus == "IGH")$cell_id
      light_cells = dplyr::filter(vdj_data, locus == "IGL" | locus == "IGK")$cell_id
      paired = intersect(heavy_cells, light_cells)
      vdj_data = dplyr::filter(vdj_data, cell_id %in% paired)
    }
    gc()
  }
  return(vdj_data)
}
