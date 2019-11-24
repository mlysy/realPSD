#' Bin a long vector and return the bin average values.
#'
#' Compute average values of a series in pre-determined bins (column-wise subsets).
#'
#' @param x Numeric vector.
#' @param bin_size Size of each bin.
#' @return A vector of length \code{floor(length(x)/bin_size) * bin_size} containing the bin averages.
#' @export
binning <- function(x, bin_size) {
  nbins <- floor(length(x) / bin_size) # number of bins
  colMeans(matrix(x[1:(nbins*bin_size)], nrow = bin_size))
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
