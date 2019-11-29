#' Simulate a Gaussian time-domain signal from its continuous-time PSD.
#'
#' @param fs Sampling frequency.
#' @param N Total number of observations.
#' @param fseq Vector of (ordered) frequencies at which the continuous-time PSD is provided.
#' @param cpsd Vector of continuous-time PSD values corresponding to `fseq`.
#' @return Simulated time series (real numeric vector).
#'
#' @details
#' @export
sim_time <- function(fs, N, fseq, cpsd) {
  # check frequency basis
  if(fseq[length(fseq)] < fs/2) {
    stop("max(fseq) must be greater than fs/2.")
  }
  if(fseq[1] > fs/N) {
    stop("min(fseq) must be less than fs/N.")
  }
  yfreq <- rep(NA_complex_, N) # storage
  # DC offset
  if(fseq[1] == 0 & cpsd[1] < Inf) {
    yfreq[1] <- sqrt(fs * cpsd[1]) * rnorm(1)
  } else yfreq[1] <- 0
  if(N == 1) return(Re(yfreq))
  # discrete-time psd: approximated as fs * cpsd
  dpsd <- fs * approx(x = fseq, y = cpsd,
                      xout = fft_fseq(fs = fs, N = N))$y
  # simulation in complex domain
  nfreq <- length(dpsd)
  yfreq[1+1:nfreq] <- sqrt(dpsd/2) * (rnorm(nfreq) + 1i * rnorm(nfreq))
  # Nyquist frequency
  if(N_even <- (N%%2) == 0) yfreq[1+nfreq] <- sqrt(2) * Re(yfreq[1+nfreq])
  # unfold
  yfreq[nfreq+(!N_even)+1:nfreq] <- Conj(yfreq[1+nfreq:1])
  # convert to time domain
  Re(fftw::IFFT(yfreq * sqrt(N)))
  ## yfreq <- sqrt(dpsd/2)
  ## check_fseq(fs, N, fseq)
  ## # sampling frequency
  ## N <- 2*N # double time series length, throw out
  ## # error handling
  ## if(fseq[length(fseq)]*2 < fs)
  ##   stop("Sampling rate exceeds the range of the spectrum.")
  ## if(1/fseq[1+ (fseq[1] == 0)] < (N-1)/fs)
  ##   message(paste0("Length of time series exceeds spectral period.\n",
  ##     "Lower frequencies taken to be minimum frequency value."))
  ## # generate frequencies
  ## m <- floor(N/2)
  ## fTs <- (1:m) * (fs/N) # Frequency range of the time series
  ## PTs <- approx(fseq, cpsd, fTs, method = "linear", cpsd[1], cpsd[1]) # Spectrum of the time series
  ## PTs <- PTs$y # extract the spectrum values
  ## # simulate a Fourier basis
  ## yFreq <- rnorm(length(fTs)) + 1i*rnorm(length(fTs))
  ## # One sqrt(2) is for unfolding, the other is for the mean of the sum of
  ## # squares of two iid Gaussians
  ## yFreq <- yFreq * sqrt(PTs)/2
  ## # Unfold (Nyquist frequency doesn't undergo aliasing and is real)
  ## n <- length(yFreq)
  ## yFreq[n] <- Re(yFreq[n])*sqrt(2)
  ## yConj <- Conj(rev(yFreq[1:(n-1)]))
  ## # DC offset
  ## dcOffset <- (fseq[1] == 0) * cpsd[1]
  ## # combine them together
  ## yFreq <- c(dcOffset, yFreq, yConj)
  ## # Recover the time series
  ## FR <- fs/N;
  ## yTime <- fftw::IFFT(yFreq*sqrt(FR)*N)
  ## yTime <- Re(yTime[1: (length(yTime)/2)])
  ## xTime <- (0: (length(yTime)-1))/fs
  ## # return
  ## return(list(xTime = xTime, yTime = yTime))
}

#--- helper functions ----------------------------------------------------------

#--- depreciated ---------------------------------------------------------------

## tsSim <- function(SF, N, xPSD, yPSD) {
##   # sampling frequency
##   N <- 2*N # double time series length, throw out
##   # error handling
##   if(xPSD[length(xPSD)]*2 < SF)
##     stop("Sampling rate exceeds the range of the spectrum.")
##   if(1/xPSD[1+ (xPSD[1] == 0)] < (N-1)/SF)
##     message(paste0("Length of time series exceeds spectral period.\n",
##                    "Lower frequencies taken to be minimum frequency value."))
##   # generate frequencies
##   m <- floor(N/2)
##   fTs <- (1:m) * (SF/N) # Frequency range of the time series
##   PTs <- approx(xPSD, yPSD, fTs, method = "linear", yPSD[1], yPSD[1]) # Spectrum of the time series
##   PTs <- PTs$y # extract the spectrum values
##   # simulate a Fourier basis
##   yFreq <- rnorm(length(fTs)) + 1i*rnorm(length(fTs))
##   # One sqrt(2) is for unfolding, the other is for the mean of the sum of
##   # squares of two iid Gaussians
##   yFreq <- yFreq * sqrt(PTs)/2
##   # Unfold (Nyquist frequency doesn't undergo aliasing and is real)
##   n <- length(yFreq)
##   yFreq[n] <- Re(yFreq[n])*sqrt(2)
##   yConj <- Conj(rev(yFreq[1:(n-1)]))
##   # DC offset
##   dcOffset <- (xPSD[1] == 0) * yPSD[1]
##   # combine them together
##   yFreq <- c(dcOffset, yFreq, yConj)
##   # Recover the time series
##   FR <- SF/N;
##   yTime <- fftw::IFFT(yFreq*sqrt(FR)*N)
##   yTime <- Re(yTime[1: (length(yTime)/2)])
##   xTime <- (0: (length(yTime)-1))/SF
##   # return
##   return(list(xTime = xTime, yTime = yTime))
## }


