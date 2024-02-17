#' Fetch all genes with specific prefix expressed in the data
#' This function fetches all genes with the specified prefix expressed in the dataset
#' @param seurat_object A Seurat object
#' @param gene_prefix The prefix of gene names to be retrieved
#' @export
#' @examples
#' KG_fetch_genes(pbmc_dataset, "HLA-")
KG_fetch_genes = function(seurat_object, gene_prefix){
  gene_list = grep(pattern = gene_prefix,
                   x = rownames(seurat_object),
                   value = T)
  return(gene_list)
}
