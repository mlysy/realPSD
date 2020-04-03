#' Calculate the bin size correction factor for the `LP` estimator.
#'
#' @param B Bin size (integer).
#'
#' @return The bin size correction factor `log(B) - digamma(B)`.
#' @export
bin_factor <- function(B) log(B) - digamma(B)
