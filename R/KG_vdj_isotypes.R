#' Annotate isotypes in V(D)J-seq data
#' 
#' This function adds a column in a data frame with V(D)J-seq data for isotypes absed on contents in the c_call column.
#' @param bcr_data A data frame containing the BCR data from V(D)J-seq.
#' @param c_call Name of the column where C calls are located.
#' @export
#' @examples
#' KG_vdj_isotypes(BCR_object)
KG_vdj_isotypes = function(bcr_data, 
                           c_call = "c_call"){
  bcr_data[["isotype"]] = NA
  
  for(i in 1:nrow(bcr_data)){
    if(bcr_data[[c_call]][i] %in% "IGHM"){
      bcr_data[["isotype"]][i] = "IgM"
    }
    else if(bcr_data[[c_call]][i] %in% c("IGHA1", "IGHA2")){
      bcr_data[["isotype"]][i] = "IgA"
    } 
    else if(bcr_data[[c_call]][i] %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4")){
      bcr_data[["isotype"]][i] = "IgG"
    } 
    else if(bcr_data[[c_call]][i] %in% "IGHE"){
      bcr_data[["isotype"]][i] = "IgE"
    } 
    else if(bcr_data[[c_call]][i] %in% "IGHD"){
      bcr_data[["isotype"]][i] = "IgD"
    } 
    else if(bcr_data[[c_call]][i] %in% "IGKC"){
      bcr_data[["isotype"]][i] = "IgK"
    } 
    else if(bcr_data[[c_call]][i] %in% c("IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7")){
      bcr_data[["isotype"]][i] = "IgL"
    } 
    else {
      bcr_data[["isotype"]][i]
    }
  }
  return(bcr_data)
}