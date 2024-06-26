% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KG_plot_to_pdf.R
\name{KG_plot_to_pdf}
\alias{KG_plot_to_pdf}
\title{Print plots trees to PDF.}
\usage{
KG_plot_to_pdf(
  plots,
  ncol = 3,
  nrow = 5,
  heading = NULL,
  padding = NULL,
  file = NULL,
  pdf_size = a4,
  portrait = TRUE
)
}
\arguments{
\item{plots}{A list of plots.}

\item{ncol}{The number of columns per page.}

\item{nrow}{The number of rows per page.}

\item{heading}{The heading to be printed on each page.}

\item{padding}{Distance between individual plots.}

\item{file}{The name of the PDF file.}

\item{pdf_size}{Size of the PDF.}

\item{portrait}{Wheter or not to use portrait layout (if FALSE, will use horizontal layout).}
}
\description{
This function prints a list of plots to a PDF file.
}
\examples{
KG_trees_to_pdf(plots, file = "plot_file")
}
