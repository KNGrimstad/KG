#' Convert isotypes to Switched and Unswitched
#'
#' This function annotates the isotype switch status (unswitched vs. switched) based on annotated isotypes in a Seurat object.
#' @param seurat_object A Seurat object
#' KG_isotype_switch(B_cell_dataset)
KG_isotype_switch = function(seurat_object){
  # Check that the isotype slot exists
  if(!("isotype" %in% names(seurat_object@meta.data))){
    stop(paste("Meta data slot isotype not found in the Seurat object.\n Please run KG_gene_to_isotype."))
  }
  isotypes = as.character(seurat_object$isotype)
  for(i in 1:ncol(seurat_object)){
    if(isotypes[i] %in% c("IgA", "IgG", "IgE")){
      isotypes[i] = "Switched"
    } else if(isotypes[i] %in% c("IgM", "IgD")){
      isotypes[i] = "Unswitched"
    } else{
      isotypes[i] = isotypes[i]
    }
  }
  seurat_object[["switch"]] = factor(isotypes,
                                     levels = c("Unswitched", "Switched"))
  return(seurat_object)
}
