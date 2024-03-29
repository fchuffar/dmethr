% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common.R
\name{et_gsea_plot}
\alias{et_gsea_plot}
\title{Draw enrichment plot}
\usage{
et_gsea_plot(expression_vector, gene_set, prefix, nperm = 1000,
  PLOT_GSEA = FALSE)
}
\arguments{
\item{expression_vector}{named numeric vector}

\item{gene_set}{character vector of gene of interest}

\item{prefix}{string used to prefix outputs}

\item{nperm}{interger number of permutation to use}

\item{PLOT_GSEA}{boolean defining if the GSEA plot needs to be computed}
}
\description{
This function enrichment plot.
}
