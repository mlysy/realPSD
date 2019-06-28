#' @title Simulates a time series from a given PSD at a given sampling frequency.
#' 
#' @param SF Sampling frequency.
#' @param N Total number of samples.
#' @param xPSD Frequency basis, possibly including DC offset (xPSD(1) = 0) if DC offset is not provided then it is set to 0.
#' @param yPSD Power spectrum up to Nyquist frequency.
#' @return A list of xTime (Time observation times) and yTime (Time series, i.e., amplitude).
#' @details ...
#' @export
tsSim <- function(SF, N, xPSD, yPSD) {
  # sampling frequency
  N <- 2*N # double time series length, throw out 
  # error handling
  if(xPSD[length(xPSD)]*2 < SF)
    stop("Sampling rate exceeds the range of the spectrum.")
  if(1/xPSD[1+ (xPSD[1] == 0)] < (N-1)/SF)
    message(paste0("Length of time series exceeds spectral period.\n",
      "Lower frequencies taken to be minimum frequency value."))
  # generate frequencies
  m <- floor(N/2)
  fTs <- (1:m) * (SF/N) # Frequency range of the time series
  PTs <- approx(xPSD, yPSD, fTs, method = "linear", yPSD[1], yPSD[1]) # Spectrum of the time series
  PTs <- PTs$y # extract the spectrum values
  # simulate a Fourier basis
  yFreq <- rnorm(length(fTs)) + 1i*rnorm(length(fTs))
  # One sqrt(2) is for unfolding, the other is for the mean of the sum of
  # squares of two iid Gaussians
  yFreq <- yFreq * sqrt(PTs)/2
  # Unfold (Nyquist frequency doesn't undergo aliasing and is real)
  n <- length(yFreq)
  yFreq[n] <- Re(yFreq[n])*sqrt(2)
  yConj <- Conj(rev(yFreq[1:(n-1)]))
  # DC offset
  dcOffset <- (xPSD[1] == 0) * yPSD[1]
  # combine them together
  yFreq <- c(dcOffset, yFreq, yConj)
  # Recover the time series
  FR <- SF/N;
  yTime = fftw::IFFT(yFreq*sqrt(FR)*N)
  yTime = Re(yTime[1: (length(yTime)/2)])
  xTime = (0: (length(yTime)-1))/SF
  # return
  return(list(xTime = xTime, yTime = yTime))
}
