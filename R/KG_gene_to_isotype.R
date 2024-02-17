#' Convert C calls to isotype
#'
#' This function annotates the major isotypes for each cell in a Seurat object based on the heavy chain c calls from the Cell Ranger output.
#' @param seurat_object A Seurat object.
#' @param c_call The column in the Seurat object where the C calls are stored.
#' KG_gene_to_isotype(B_cell_dataset, "c_call")
KG_gene_to_isotype = function(seurat_object, c_call){

  # Check that the specified column (c_call) exists in the Seurat object
  if (!(c_call %in% names(seurat_object@meta.data))) {
    stop(paste("Meta data slot", c_call, "not found in the Seurat object"))
  }

  # Create new column for modified isotype values
  c_calls_modified = character(length = ncol(seurat_object))

  # Loop
  for(i in 1:ncol(seurat_object)){
    gene_value = seurat_object[[c_call]][i,]
    # Check for the conditions and replace strings accordingly
    if(gene_value %in% c("IGHA1", "IGHA2")){
      c_calls_modified[i] = "IgA"
    } else if(gene_value %in% "IGHM"){
      c_calls_modified[i] = "IgM"
    } else if(gene_value %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4")){
      c_calls_modified[i] = "IgG"
    } else if(gene_value %in% "IGHD"){
      c_calls_modified[i] = "IgD"
    } else if(gene_value %in% "IGHE"){
      c_calls_modified[i] = "IgE"
    } else{
      # If none of the conditions are match, keep original string
      c_calls_modified[i] = gene_value
    }
  }
  # Add new column to the Seurat object with specialized name
  seurat_object[["isotype"]] = factor(c_calls_modified,
                                      levels = c("IgM", "IgD", "IgG", "IgA", "IgE"))
  return(seurat_object)
}
