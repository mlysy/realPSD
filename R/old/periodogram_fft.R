#' @title Calculates the one-sided periodogram from given FFT.
#' 
#' @param yFFT FFT of the original time series.
#' @param SF_s Sampling frequecy (Hz).
#' @param T_s Total time (seconds).
#' @return List of xFreq: One-sided periodogram frequencies and yFreq: One-sidede periodogram ordinates.
#' @export
periodogram_fft <- function(yFFT, SF_s, T_s) {
  # fft(yFFT)
  N <- length(yFFT) - (length(yFFT) %% 2) # use even PSD's only
  FR <- SF_s/N
  yFFT <- yFFT[1:N]
  yFreq <- yFFT/(sqrt(FR)*N)
  # Remove DC component, select first frequency up to the Nyquist frequency.
  # Drop negative frequencies. For a real signal, negative frequencies are redundant.
  yFreq <- yFreq[2:(N/2+1)]
  xFreq <- seq(from = 1/T_s, to = SF_s/2, length.out = length(yFreq))

  yFreq[length(yFreq)] <- Re(yFreq[length(yFreq)]) # Nyquist freq should be real for real signal
  yFreq <- abs(yFreq)^2 # Take the squared magnitude at each frequency
  # Make the PSD single-sided by doubling everything except the Nyquist frequency.
  # This corrects for the loss in amplitude caused by dropping negative frequencies
  # the Nyquist frequency at N/2 doesn't undergo folding, and therefore is NOT doubled.
  yFreq[1: (N/2-1)] <- 2 * yFreq[1: (N/2-1)]

  return(list(xFreq = xFreq, yFreq = yFreq))
}
