#' Discrete Fourier transform of a sine wave.
#'
#' @param freq Frequency vector.
#' @param fs Sampling frequency.
#' @param N Number of observations (see 'Details').
#' @param A Amplitude of the sine wave (see 'Details').
#' @param xi Frequency (i.e., inverse-period) of the sine wave (see 'Details').
#' @param phi Phase of the sine wave (see 'Details').
#' @return A vector of the same length at `freq` containing the continuously-extended discrete Fourier transform of the sine wave evaluated at time points `0:(N-1)/fs`.
#'
#' @details For a sine wave defined as
#' ```
#' g(t) = A * sin(2*pi*xi * t + phi),
#' ```
#' computes the (continuously-extended) discrete Fourier transform
#' ```
#' G(f) = sum_{n=0}^{N-1} g(n/fs) * exp(2*pi*i/fs * f)
#' ```
#' at the frequencies in `freq`.  When `freq = 0:(N-1) * fs/N`, `G(freq)` is the FFT of `g(0:(N-1)/fs)`.
#' @export
sin_fft <- function(freq, fs, N, A, xi, phi) {
  ## dT <- 1/fs
  ifac <- 2*pi*1i/fs
  pterm <- (exp(ifac * N * (xi-freq)) - 1) / (exp(ifac * (xi-freq)) - 1)
  mterm <- (exp(-ifac * N * (xi+freq)) - 1) / (exp(-ifac * (xi+freq)) - 1)
  A/(2*1i) * (exp(phi*1i) * pterm - exp(-phi*1i) * mterm)
}
