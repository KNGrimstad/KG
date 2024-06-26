% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KG_cluster_multi.R
\name{KG_cluster_multi}
\alias{KG_cluster_multi}
\title{Evaluate cluster resolutions with clustree}
\usage{
KG_cluster_multi(
  seurat_object,
  max.res = 1,
  min.res = 0.1,
  increment = 0.1,
  plot = TRUE,
  assay = "RNA",
  dims = 1:35,
  reduction = "pca",
  verbose = TRUE
)
}
\arguments{
\item{seurat_object}{A Seurat object.}

\item{max.res}{The maximum resolution for clustering.}

\item{min.res}{The minimum resolution for clustering.}

\item{increment}{A value determining by how much each increment in resolution should be.}

\item{plot}{If TRUE, returns a tree plot of resulting clusters for each resolution.}

\item{assay}{Name of the Seurat assay to use, "RNA" by default.}

\item{dims}{The number of dimensions to use for FindNeighbors.}

\item{reduction}{Name of the dimensionality reduction method to use, "pca" by default.}
}
\description{
This function runs the function FindClusters from the Seurat package, for a range of resolution values, plots the clusters for each resolution, and prompts the user to select an appropriate resolution, and returns the clustered Seurat object for the selected resolution.
}
\examples{
KG_cluster_multi(B_cell_dataset)
}
