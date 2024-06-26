% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KG_isotype_dimplot.R
\name{KG_isotype_dimplot}
\alias{KG_isotype_dimplot}
\title{Make a nice dimplot}
\usage{
KG_isotype_dimplot(
  seurat_object,
  cells = NULL,
  reduction = "umap",
  dims = c(1, 2),
  group.by = NULL,
  label = TRUE,
  label.size = 7,
  pt.size = 2,
  stroke = 0.1,
  legend = FALSE,
  axes = TRUE,
  heading = TRUE,
  isotype = "isotype",
  combine = TRUE,
  ncol = 2
)
}
\arguments{
\item{seurat_object}{A Seurat object.}

\item{cells}{List of cells to plot. If NULL, all cells are plotted.}

\item{reduction}{The dimensionality reduction to use.}

\item{dims}{Dimensions to plot. By default, plots the first two dimensions.}

\item{group.by}{Group identifier to use for shading.}

\item{label}{Whether to plot group identifier labels.}

\item{label.size}{Size of the group identifier labels.}

\item{pt.size}{Size of data points.}

\item{stroke}{Size of border to data points. NOTE: when exportet as PDF, borders can alter the appearance of the plot. To avoid this, export externally as PNG file.}

\item{legend}{Whether or not to plot the legend/key to the group identifier.}

\item{axes}{Logical; whether or not to plot the axes.}

\item{heading}{Logical; whether or not to include a heading (with the name of the isotype plotted).}

\item{isotype}{Slow in Seurat object where the isotype information is stored.}

\item{combine}{Logical; whether or not to print all plots on a single page.}

\item{ncol}{Number of columns to organize the plots in when printed.}
}
\description{
This function constructs a nice, clean dimplot from dimensionally reduced data points in a Seurat object.
}
\examples{
KG_isotype_dimplot(B_cell_dataset)
}
