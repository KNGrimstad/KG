% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KG_percent_cells.R
\name{KG_percent_cells}
\alias{KG_percent_cells}
\title{Plot percentage of cells for each cluster.}
\usage{
KG_percent_cells(
  seurat_object,
  group.by = NULL,
  plot = "pie",
  cols = NULL,
  label = T,
  label_size = 3,
  label_col = "black",
  legend = T
)
}
\arguments{
\item{seurat_object}{A Seurat object}

\item{group.by}{The identity column to use from the Seurat object.}

\item{plot}{How to plot the results. Options include "table", "pie", "bar".}

\item{cols}{Vector of colors to use.}

\item{label}{Logical; whether or not to add the percentages as text albels in the plot. Irrelevant for tables.}

\item{label_size}{The size of labels.}

\item{label_col}{Color for the labels.}

\item{legend}{Logical; whether or not to include legend.}
}
\description{
This function plots the percentage of cells, out of the entire dataset, that each cluster contains.
}
\examples{
KG_percent_cells(B_cell_dataset)
}
