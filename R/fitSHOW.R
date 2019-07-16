# fit simulated datasets
#' @param f Vector of frequencies, usually from f0 - f0/sqrt(2) to f0 + f0/sqrt(2).
#' @param rfreq Vector of exponential random variables Exp(1) with the same length as f.
#' @param fs Sampling frequency, Hz.
#' @param f0 Resonance frequency, Hz.
#' @param Q Quality factor.
#' @param k Cantilever stiffness, N/m.
#' @param Kb Boltzmann's constant.
#' @param T Temperature, Kelvin.
#' @param Aw White noise psd.
#' @param binSize Integer number, bin size.
#' @param method Fitting method, i.e. LP_nlp, MLE_nlp, NLS_nlp.
#' @export
fitSHOW <- function(f, rfreq, fs, f0, Q, k, T, Aw, 
    binSize = 100, method = c("LP_nlp", "MLE_nlp", "NLS_nlp")) {
  Kb <- 1.381e-23           # Boltzmann's constant
  sig2 <- Kb*T/(k*pi*f0*Q) # variance
  Rw <- Aw/sig2 # re-parameterization
  phi <- c(f0, f0*Q, Rw) # parameter vector for SHOW model
  # psd values at each frequency point of f with given Q
  psd <- psdSHO(f, f0, Q, k, Kb, T, unit_conversion = TRUE) + Aw
  # generate the periodogram values
  Y <- rfreq * fs * psd 
  # ---------- binning ----------
  # binSize <- 100
  fbar <- binning(f, binSize = binSize)
  Ybar <- binning(Y, binSize = binSize)
  Zbar <- log(Ybar)
  # ---------- fitting ----------
  if (method == "LP_nlp") {
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = method,
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(phi)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
  } else if (method == "NLS_nlp") {
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = method,
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar)),
                           parameters = list(phi = matrix(phi)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
  } else if (method == "MLE_nlp") {
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = method,
                                       f = matrix(f),
                                       Y = matrix(Y)),
                           parameters = list(phi = matrix(phi)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
  } else {
    stop("method should be chosen from LP_nlp, NLS_nlp and MLE_nlp.")
  }
  opt <- optim(tmod$par, fn = tmod$fn, gr = tmod$gr, control = list(maxit = 1000))
  phi_hat <- opt$par # extract the fitted parameters
  param <- rep(NA, 3) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2]/phi_hat[1] # Q_hat
  param[3] <- phi_hat[3] / Aw * Kb * T / (pi * phi_hat[2]) # k_hat
  names(param) <- c("f0", "Q", "k")
  return(param)
}
