#' @title binning
#' @description
#' Compute average values of a series in pre-determined bins (column-wise subsets).
#' The bin size can be determined either directly or by specifying the number of bins.
#' @param x numeric \code{vector} to process
#' @param binSize length of each bin
#' @return a \code{matrix} or \code{vector} with average values per bin
#' @export
binning <- function(x, binSize) {
  # number of bins 
  nbins <- ceiling(length(x) / binSize)
  # number of positions left
  nleft <- nbins*binSize - length(x) 
  # fill the remaining positions with 0s
  tmp <- c(x, rep(0, nleft)) 
  # bin the vector and calculate the average values
  mat <- matrix(tmp, nrow = binSize) # each column is a bin 
  xbar <- colMeans(mat)
  # correct the last bin average
  xbar[nbins] <- sum(x[(binSize*(nbins-1)+1):length(x)])/(length(x) - binSize*(nbins-1))
  # return 
  return(xbar)
} 
