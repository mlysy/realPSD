#' @title binning a long vector and return the bin average values
#' @description
#' Compute average values of a series in pre-determined bins (column-wise subsets).
#' The bin size can be determined either directly or by specifying the number of bins.
#' @param x numeric \code{vector} to process
#' @param bin_size length of each bin
#' @return a \code{matrix} or \code{vector} with average values per bin
#' @export
binning <- function(x, bin_size) {
  # number of bins 
  nbins <- ceiling(length(x) / bin_size)
  # number of positions left
  nleft <- nbins*bin_size - length(x) 
  # fill the remaining positions with 0s
  tmp <- c(x, rep(0, nleft)) 
  # bin the vector and calculate the average values
  mat <- matrix(tmp, nrow = bin_size) # each column is a bin 
  xbar <- colMeans(mat)
  # correct the last bin average
  xbar[nbins] <- sum(x[(bin_size*(nbins-1)+1):length(x)])/(length(x) - bin_size*(nbins-1))
  # return 
  return(xbar)
} 
