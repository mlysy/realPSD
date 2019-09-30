#' @title Fit simulated datasets to get fitted parameters (for internal use)
#' 
#' @param fseq Sequence of frequencies, usually from f0 - f0/sqrt(2) to f0 + f0/sqrt(2).
#' @param sim_cnorm Vector of complex normal random variables 
#' @param f0 Resonance frequency, Hz.
#' @param fs Sampling frequency, Hz.
#' @param Q Quality factor.
#' @param k Cantilever stiffness, N/m.
#' @param Temp Temperature, Kelvin.
#' @param add_white_noise TRUE/FALSE indicator
#' @param Aw White noise psd.
#' @param Nfreq Total number of frequency points in the whole range.
#' @param bin_size Integer number, bin size.
#' @param method Fitting method, i.e. lp, mle, nls.
#' @param unit_conversion Logical, if TRUE, use fm2/Hz instead of m2/Hz
fitSHOWsine <- function(fseq, sim_cnorm, f0, fs, Q, k, Temp, Aw, Nfreq,
                    add_white_noise = TRUE,
                    bin_size = 100, method = c("lp", "mle", "nls"),
                    remove_noise = TRUE,
                    unit_conversion = FALSE) {
  # ---------- setup -----------
  method <- match.arg(method)
  Kb <- 1.381e-23           # Boltzmann's constant
  if(unit_conversion) {
    sig2 <- Kb*Temp/(k*pi*f0*Q) * 1e30 # variance, unit: fm2/Hz 
  } else {
    sig2 <- Kb*Temp/(k*pi*f0*Q) # variance, unit: m2/Hz
  }
  Rw <- Aw/sig2 # re-parameterization, note we input Aw with unit fm2/Hz
  phi <- c(f0, Q, Rw) # parameter vector for SHOW model
  # phi <- c(f0 + rnorm(1, 0, sqrt(f0)/10), 
  #   f0*Q + rnorm(1, 0, sqrt(f0*Q)/10), 
  #   Rw + rnorm(1,0, Rw/10)) 
  # psd values at each frequency point of f with given Q
  if(add_white_noise) {
    psd <- psdSHO(fseq, f0, Q, k, Temp, unit_conversion) + Aw
  } else {
    psd <- psdSHO(fseq, f0, Q, k, Temp, unit_conversion)
  }
  # generate the periodogram values
  sin_fft <- fft_sin(fseq, Nfreq, f0, Q, fs, unit_conversion)
  Y <- sim_cnorm * sqrt(psd * fs)
  Y <- (Y + sin_fft) * Conj(Y + sin_fft)
  Y <- Re(Y)
  # remove sine wave noise 
  if(remove_noise) {
    # convert Y to standard unit (otherwise the NLS optim would fail)
    if(unit_conversion) Y <- Y/1e30
    # preliminary estimation
    param_pre <- fitSHOW_TMB(fseq, Y, bin_size, method, phi, Temp, Kb)
    # if optim above returns NA, then skip this whole estimation
    if(any(is.na(param_pre))) return(rep(NA,4))
    # remove sine wave noise
    freq_range <- c(f0-f0/sqrt(2), f0+f0/sqrt(2))
    f0_hat <- param_pre["f0_hat"]
    Q_hat <- param_pre["Q_hat"]
    k_hat <- param_pre["k_hat"]
    Aw_hat <- param_pre["Aw_hat"]
    if(unit_conversion) Y <- Y * 1e30
    Y <- psd_denoise(fseq, psd_noise = Y, 
      Q_hat, f0_hat, k_hat, Temp, unit_conversion, Aw_hat, freq_range)
  }
  # convert Y to standard unit (otherwise the NLS optim would fail)
  if(unit_conversion) Y <- Y/1e30
  # param estimation
  param <- fitSHOW_TMB(fseq, Y, bin_size, method, phi, Temp, Kb)
  return(param)
}
