#' Calculate the one-sided periodogram of a time series.
#'
#' @param Xt Vector of time series observations.
#' @param fs Sampling frequecy (Hz).
#' @return Data frame with columns `freq` and `Ypsd`.
#'
#' @details For `length(Xt) = N = 2*K`, we have
#' ```
#' freq = 1:(K-1)/N * fs,     Ypsd = 1/(N*fs) * abs(fft(Xt)[2:K])^2.
#' ```
#' @export
periodogram <- function(Xt, fs) {
  # fft(Xt)
  N <- length(Xt) - (length(Xt) %% 2) # use even PSD's only
  K <- N/2
  ## FR <- fs/N
  Xt <- Xt[1:N]
  ## Ypsd <- fftw::FFT(Xt)/(sqrt(FR)*N)
  Xf <- fftw::FFT(Xt)/sqrt(fs * N)
  # Remove DC component, Nyquist frequency, and negative frequencies.
  Xf <- Xf[2:K]
  freq <- 1:(K-1) * fs/N
  Ypsd <- abs(Xf)^2 # periodogram is square magnitude of FFT
  data.frame(freq = freq, Ypsd = Ypsd)
}
