% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binning.R
\name{binning}
\alias{binning}
\title{Bin a long vector and return the bin average values.}
\usage{
binning(x, bin_size, bin_type = c("mean", "median"))
}
\arguments{
\item{x}{Numeric vector.}

\item{bin_size}{Size of each bin.}

\item{bin_type}{Type of binning (mean or median).}
}
\value{
A vector of length \code{floor(length(x)/bin_size) * bin_size} containing the bin means or medians.
}
\description{
Compute average values of a series in pre-determined bins (column-wise subsets).
}
