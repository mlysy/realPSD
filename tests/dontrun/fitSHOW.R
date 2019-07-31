#' @title Fit simulated datasets to get fitted parameters (for internal use)
#' 
#' @param fseq Sequence of frequencies, usually from f0 - f0/sqrt(2) to f0 + f0/sqrt(2).
#' @param sim_exp Vector of exponential random variables Exp(1) with the same length as fseq.
#' @param f0 Resonance frequency, Hz.
#' @param fs Sampling frequency, Hz.
#' @param Q Quality factor.
#' @param k Cantilever stiffness, N/m.
#' @param Temp Temperature, Kelvin.
#' @param Aw White noise psd.
#' @param bin_size Integer number, bin size.
#' @param method Fitting method, i.e. lp, mle, nls.
fitSHOW <- function(fseq, sim_exp, f0, fs, Q, k, Temp, Aw,
                    bin_size = 100, method = c("lp", "mle", "nls")) {
  # ---------- setup -----------
  method <- match.arg(method)
  Kb <- 1.381e-23           # Boltzmann's constant
  sig2 <- Kb*Temp/(k*pi*f0*Q) # variance, unit: m2/Hz
  Rw <- Aw/sig2 # re-parameterization, note we input Aw with unit fm2/Hz
  phi <- c(f0, f0*Q, Rw) # parameter vector for SHOW model
  # psd values at each frequency point of f with given Q
  psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion = FALSE) + Aw
  # generate the periodogram values
  Y <- sim_exp * psd * fs
  # ---------- binning ----------
  # bin_size <- 100
  fbar <- binning(fseq, bin_size = bin_size)
  Ybar <- binning(Y, bin_size = bin_size)
  Zbar <- log(Ybar)
  # ---------- fitting ----------
  if (method == "lp") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_nlp",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gz <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_zeta",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
      zeta <- gz$fn(phi)
      exp(zeta)
    }
  } else if (method == "nls") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_nlp",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_tau",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
      gt$fn(phi)
    }
  } else if (method == "mle") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_nlp",
                                      f = matrix(fseq),
                                      Y = matrix(Y)),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_tau",
                                      f = matrix(fseq),
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
  phi_hat <- opt$par # extract the fitted phi = c(f0_hat, gamma_hat, Rw_hat)
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2, unit: fm^2/Hz which should be the same as Aw
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2]/phi_hat[1] # Q_hat
  param[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
  param[4] <- phi_hat[3] * tau_hat # Aw_hat
  names(param) <- c("f0_hat", "Q_hat", "k_hat", "Aw_hat")
  return(param)
}
