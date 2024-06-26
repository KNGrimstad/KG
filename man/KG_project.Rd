% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KG_project.R
\name{KG_project}
\alias{KG_project}
\title{Project data}
\usage{
KG_project(
  seurat_object,
  reference,
  reduction = "pca",
  assay = "RNA",
  reference.assay = "RNA",
  dims = 1:35,
  min.dist = 0.3,
  spread = 1,
  reduction.model = "umap",
  reference.reduction = "pca",
  reference.idents = NULL
)
}
\arguments{
\item{seurat_object}{A Seurat object.}

\item{reference}{A Seurat object from which cluster identities and UMAP coordinates/embeddings should be pulled.}

\item{reduction}{Name of the dimensional reduction method to use as input, "pca" by default.}

\item{assay}{Name of the Seurat object assay to use.}

\item{reference.assay}{Name of the Seurat reference object assay to use.}

\item{dims}{The range of dimensions to use from linear dimensional reduction method.}

\item{min.dist}{Controls how tightly embeddings are allowed to compress points together. Larger values result in more even distribution. Smaller values optimizes more accurately the local structure.}

\item{spread}{Effective scale of embedded points. Determines, along with min.dist, how clumped embedded points are.}

\item{reduction.model}{Name of dimensional reduction model to use for plotting, "umap" by default.}

\item{reference.reduction}{Name of dimensional reduction method to use as input from reference, "pca" by default.}

\item{reference.idents}{The identity variable to use from the reference for the projection. By default, uses the default identity of the reference object.}
}
\description{
This function projects a Seurat object onto the clusters and UMAP coordinates from a reference.
}
\examples{
KG_project(B_cell_dataset, B_cell_reference_dataset)
}
