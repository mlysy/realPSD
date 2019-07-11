#' @title binning
#' @description
#' Compute average values of a series in pre-determined bins (column-wise subsets).
#' The bin size can be determined either directly or by specifying the number of bins.
#' @param X numeric \code{data.frame}, \code{matrix} or \code{vector} to process
#' @param nbins number of bins
#' @return a \code{matrix} or \code{vector} with average values per bin
#' @export
binning <- function(X, nbins) {
    tapply(X, cut(X, nbins), mean) 
    # if (is.data.frame(X)) 
    #     X <- as.matrix(X)
    # if (!missing(nbins) & !missing(binSize)) 
    #     stop("EITHER 'bins' OR 'binSize' must be specified")
    # if (missing(nbins) & missing(binSize)) 
    #     return(X)
    
    # if (is.matrix(X)) 
    #     p1 <- ncol(X) else p1 <- length(X)
    
    # if (missing(nbins) & !missing(binSize)) {
    #     b <- findInterval(1:p1, seq(1, p1, binSize))
    # } else {
    #     b <- findInterval(1:p1, seq(1, p1, length.out = nbins + 1), rightmost.closed = T)
    # }
    
    # p2 <- max(b)
    
    # if (is.matrix(X)) {
    #     output <- matrix(0, nrow(X), p2)
    #     for (i in seq_len(p2)) {
    #         output[, i] <- rowMeans(X[, b == i, drop = F])
    #     }
    #     colnames(output) <- colnames(X)[ceiling(tapply(b, b, function(x) mean(which(b == x[1]), na.rm = T)))]  # find colnames
    #     rownames(output) <- rownames(X)
    # } else {
    #     output <- tapply(X, b, mean)
    #     names(output) <- names(X)[ceiling(tapply(b, b, function(x) mean(which(b == x[1]), na.rm = T)))]
    # }
    
    # return(output)
} 
