#' @title fit simulated datasets to get fitted parameters
#' @param fseq Sequence of frequencies, usually from f0 - f0/sqrt(2) to f0 + f0/sqrt(2).
#' @param sim_exp Vector of exponential random variables Exp(1) with the same length as fseq.
#' @param fs Sampling frequency, Hz.
#' @param f0 Resonance frequency, Hz.
#' @param Q Quality factor.
#' @param k Cantilever stiffness, N/m.
#' @param Temp Temperature, Kelvin.
#' @param Aw White noise psd.
#' @param bin_size Integer number, bin size.
#' @param method Fitting method, i.e. lp, mle, nls.
#' @export
fitSHOW <- function(fseq, sim_exp, fs, f0, Q, k, Temp, Aw,
                    bin_size = 100, method = c("lp", "mle", "nls")) {
  method <- match.arg(method)
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
  if (method == "lp") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_nlp",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar)),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gz <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "LP_zeta",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(0, 3, 1)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
      zeta <- gz$fn(phi)
      exp(zeta)
    }
  } else if (method == "nls") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_nlp",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar)),
                           parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "NLS_tau",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar)),
                           parameters = list(phi = matrix(0, 3, 1)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
      gt$fn(phi)
    }
  } else if (method == "mle") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_nlp",
                                      f = matrix(f),
                                      Y = matrix(Y)),
                           parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "MLE_tau",
                                       f = matrix(f),
                                       Y = matrix(Y)),
                           parameters = list(phi = matrix(0, 3, 1)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
      gt$fn(phi)
    }
  } else {
    stop("method should be chosen from lp, nls and mle.")
  }
  opt <- optim(phi, fn = obj$fn, gr = obj$gr,
               control = list(maxit = 1000))
  phi_hat <- opt$par # extract the fitted parameters
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2]/phi_hat[1] # Q_hat
  param[3] <- phi_hat[3] / Aw * Kb * T / (pi * phi_hat[2]) # k_hat
  # param[4] = ???
  names(param) <- c("f0", "Q", "k", "Aw")
  return(param)
}
