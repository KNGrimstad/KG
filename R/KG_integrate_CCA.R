#' Integrate data with Seurat (CCA)
#'
#' A function for streamlining integration of Seurat objects with Seurat's CCA model.
#' @param seurat_objects A list of Seurat objects to integrate.
#' @param assay Name of the assay to use.
#' @param nIntegration_features The number of integration features to use.
#' @param reduction Name of dimensional reduction method to use in integration
#' @param dims The range of principal components to use.
#' @param scale_data Logical, whether or not to scale the data following integration.
#' @param new_name The name to give to the integrated assay.
#' @export
#' @examples
#' KG_integrate_CCA(B_cell_datasets)
KG_integrate_CCA = function(seurat_objects,
                            assay = "SCT",
                            nIntegration_features = 10000,
                            reduction = "cca",
                            dims = 1:35,
                            scale_data = TRUE,
                            new.name = "Integrated"){

  # Find most variable features across samples to integrate
  cat(paste("Finding",
            nIntegration_features,
            "integration features\n",
            sep = " "))
  integ_features = SelectIntegrationFeatures(seurat_objects,
                                             nfeatures = nIntegration_features)
  gc()

  # Prep data for integration
  cat("Preparing data for integration\n")
  seurat_int = PrepSCTIntegration(seurat_objects,
                                  anchor.features = integ_features)
  gc()
  # Identify integration anchors
  cat("Finding integration anchors\n")
  anchors = FindIntegrationAnchors(seurat_int,
                                   dims = dims,
                                   reduction = reduction,
                                   anchor.features = integ_features,
                                   normalization.method = "SCT")
  gc()
  # Integrate the datasets
  cat("Integrating datasets\n")
  seurat_int = IntegrateData(anchors,
                             dims = dims, new.assay.name = new.name, normalization.method = "SCT")
  gc()
  if (scale_data == TRUE){
    message("Scaling data")
    seurat_int = ScaleData(seurat_int, features = rownames(seurat_int))
  }
  gc()
  # Run PCA
  cat("Running principal component analysis\n")
  seurat_int = RunPCA(seurat_int,
                      verbose = F)

  return(seurat_int)
}
