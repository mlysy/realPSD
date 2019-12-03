#' Bin a long vector and return the bin average values.
#'
#' Compute average values of a series in pre-determined bins (column-wise subsets).
#'
#' @param x Numeric vector.
#' @param bin_size Size of each bin.
#' @param bin_type Type of binning (mean or median).
#' @return A vector of length \code{floor(length(x)/bin_size) * bin_size} containing the bin means or medians.
#' @export
binning <- function(x, bin_size, bin_type = c("mean", "median")) {
  bin_type <- match.arg(bin_type)
  nbins <- floor(length(x) / bin_size) # number of bins
  xmat <- matrix(x[1:(nbins*bin_size)], nrow = bin_size)
  if(bin_type == "mean") {
    out <- colMeans(xmat)
  } else if (bin_type == "median") {
    out <- apply(xmat, 2, median)
  }
  out
  ## # number of positions left
  ## nleft <- nbins*bin_size - length(x)
  ## # fill the remaining positions with 0s
  ## tmp <- c(x, rep(0, nleft))
  ## # bin the vector and calculate the average values
  ## mat <- matrix(tmp, nrow = bin_size) # each column is a bin
  ## xbar <- colMeans(mat)
  ## # correct the last bin average
  ## xbar[nbins] <- sum(x[(bin_size*(nbins-1)+1):length(x)])/(length(x) - bin_size*(nbins-1))
  ## # return
  ## return(xbar)
}
