#' Calculate percentage expression of Ig isotypes
#'
#' This function calculates the percentage expression of Ig isotypes in each cluster in a Seurat object
#' @param seurat_object A Seurat object
#' @param MGA_only Whether or not to only include IgM, IgG, and IgG isotypes
#' @param return_table Whether the output should include a data frame of the percentages
#' KG_percent_isotypes(B_cell_dataset)
KG_percent_isotypes = function(seurat_object,
                               MGA_only = TRUE,
                               return_table = FALSE){
  # Initialize an empty data frame to store results
  results_df = data.frame()

  # Calculate percentages
  tab = table(Idents(seurat_object), seurat_object$isotype)
  if(MGA_only == TRUE){
    tab = tab[, c(1, 3, 4)]
    for(i in seq_len(nrow(tab))){
      percentages = c(
        IgM = sum(tab[i, 1]) / sum(tab[i,]) * 100,
        IgG = sum(tab[i, 2]) / sum(tab[i,]) * 100,
        IgA = sum(tab[i, 3]) / sum(tab[i,]) * 100
      )

      # Print cluster information
      if(return_table == FALSE){
        cat("Cluster", rownames(tab)[i], "\n")
        cat(paste(names(percentages), ": ", round(percentages, 1), "%", sep = "", collapse = ", "), "\n")
        cat("----------\n")
      }

      # Add to results data frame
      results_df = rbind(results_df, cbind(Cluster = rownames(tab)[i], round(t(as.data.frame(percentages)), 1)))
    }

  } else{
    for(i in seq_len(nrow(tab))){
      percentages = c(
        IgM = sum(tab[i, 1]) / sum(tab[i,]) * 100,
        IgD = sum(tab[i, 2]) / sum(tab[i,]) * 100,
        IgG = sum(tab[i, 3]) / sum(tab[i,]) * 100,
        IgA = sum(tab[i, 4]) / sum(tab[i,]) * 100,
        IgE = sum(tab[i, 5]) / sum(tab[i,]) * 100
      )

      # Print cluster information
      if(return_table == FALSE){
        cat("Cluster", rownames(tab)[i], "\n")
        cat(paste(names(percentages), ": ", round(percentages, 1), "%", sep = "", collapse = ", "), "\n")
        cat("----------\n")
      }

      # Add results to results data frame
      results_df = rbind(results_df, cbind(Cluster = rownames(tab)[i], round(t(as.data.frame(percentages)), 1)))
    }
  }
  if(return_table == TRUE){
    rownames(results_df) = NULL
    return(results_df)
  }
}
