#' Calculate the tail probability of Fisher's G statistic.
#'
#' @param a Probability threshold (scalar between 0 and 1).
#' @param q Number of variables from which Fisher's G is computed (see 'Details').
#' @param log_sort If `TRUE`, reorders the terms of the alternating series to minimize numerical overflow.
#' @return Tail probability of Fisher's G statistic: `Pr(G > a)`.  Attempts to fail with `NA` when numerical underflow (`prob < 0`) is detected.
#'
#' @details Let `U = U_0 < ... < U_q`, where `U_1, ... ,U_{q-1}` are the *order statistics* of `n-1` iid Uniforms, `U_0 = 0`, and `U_q = 1`.  Then Fisher's G statistic is defined as `G = max(diff(U))`.
#'
#' It turns out that if `Yf` are `n` periodogram ordinates corresponding to white noise, then `G = max(Yf)/sum(Yf)`, such that Fisher's G statistic can be used to test for hidden periodic components in `Yf`.
#' @export
fisherG <- function(a, q, log_sort = TRUE) {
  # Maximum number of values to calculate
  maxq <- min(q, floor(1/a))
  maxq2 <- ceiling(maxq/2)
  # Calculate the terms on the log scale
  logt <- rep(0, maxq2*2)
  ind <- c(1:maxq)
  logt[ind] <- lgamma(q+1) - lgamma(q-ind+1) - lgamma(ind+1)
  logt[ind] <-  logt[ind] + (q-1) * log(1-ind*a)
  if(maxq < 2*maxq2) logt[2*maxq2] <- -Inf
  # Group the +ve and -ve terms by two to calculate the sum
  logt <- matrix(logt, 2, maxq2)
  if(log_sort && ncol(logt) > 1) {
    logt <- t(apply(logt, 1, sort)) # sort each row of logt in ascending order
  }
  maxlt <- pmax(logt[1,], logt[2,]) # apply(logt, 2, max)
  tmp <- matrix(maxlt, nrow = 2, ncol = maxq2, byrow = TRUE)
  logt <- exp(logt - tmp)
  logt <- logt[1,] - logt[2,]
  # If both terms are minuscule might get Inf-Inf
  logt[is.na(logt)] <- 0
  # Add the remaining terms
  sgt <- sign(logt)
  logt <- log(abs(logt)) + maxlt
  mx <- max(logt)
  prob <- sum(sgt * exp(logt - mx))
  if(prob <= 0) return(NA)
  prob <- exp(log(prob) + mx)
  return(prob)
}
