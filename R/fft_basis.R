#' Calculate the frequency basis of the discrete Fourier transform.
#'
#' @param fs Sampling frequency.
#' @param N Number of time domain observations.
#' @param frng Optional frequency range (vector of length 2) on which to restrict the basis.
#'
#' @return The one-sided frequency basis `freq = 1:floor(N/2)/N * fs`, optionally restricted to `frng`.
#'
#' @export
fft_basis <- function(fs, N, frng) {
  df <- fs/N # frequency discretization
  if(missing(frng)) {
    irng <- c(1, floor(N/2))
    ## frng <- c(1, floor(N/2))/N * fs
  } else {
    irng <- c(ceiling(frng[1]/df), floor(frng[2]/df))
  }
  # frequency range restricted to full basis
  irng[1]:irng[2] * df
  ## seq(from = ceiling(frng[1]/df)*df, to = floor(frng[2]/df)*df, by = df)
}
