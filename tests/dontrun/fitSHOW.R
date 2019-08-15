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
  # phi <- c(f0 + rnorm(1, 0, sqrt(f0)/10), 
  #   f0*Q + rnorm(1, 0, sqrt(f0*Q)/10), 
  #   Rw + rnorm(1,0,Rw/10)) 
  # psd values at each frequency point of f with given Q
  psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion = FALSE) + Aw
  # generate the periodogram values
  Y <- sim_exp * psd * fs
  # ---------- binning ----------
  # bin_size <- 100
  fbar <- binning(fseq, bin_size = bin_size)
  Ybar <- binning(Y, bin_size = bin_size)
  Zbar <- log(Ybar)
  # bias_correction for LP method
  bias <- digamma(bin_size) - log(bin_size)
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
      # correct the bias
      exp(zeta - bias)
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
                                      Y = matrix(Y),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_tau",
                                      f = matrix(fseq),
                                      Y = matrix(Y),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
      gt$fn(phi)
    }
  } else {
    stop("method should be chosen from lp, nls and mle.")
  }
  # ---------- optimization -----------
  opt <- optim(phi, fn = obj$fn, gr = obj$gr,
            method = "BFGS",
            control = list(maxit = 2000))
  phi_hat <- opt$par
  # check convergence 
  if(opt$convergence != 0) 
    warning(paste0(method, " didn't converge!"))
  # if(method == "mle" || method == "lp"){
  #   opt <- optim(phi, fn = obj$fn, gr = obj$gr,
  #             method = "BFGS",
  #             control = list(maxit = 2000))
  #   phi_hat <- opt$par
  # } else {
  #   names(phi) <- c("f0", "gamma", "Rw")
  #   # phi[2] <- f0 * Q # set the initial gamma (i.e. Q) to be the true value
  #   # optimize gamma (i.e., Q with f0 fixed)
  #   opt1 <- nlsr::nlfb(start = phi,
  #                     resfn = obj$fn,
  #                     lower = 0,
  #                     maskdix = c(1,0,1), # indices of parameters to be fixed
  #                     weights = rep(1,3),
  #                     control = list(jemax = 2000))
  #   # optimize f0 and Q
  #   opt2 <- nlsr::nlfb(start = opt1$coefficients,
  #                      resfn = obj$fn,
  #                      lower = 0,
  #                      maskdix = c(0,0,1),
  #                      weights = rep(1,3),
  #                      control = list(jemax = 2000))
  #   # optimize all three parameters
  #   opt3 <- nlsr::nlfb(start = opt2$coefficients,
  #                      resfn = obj$fn,
  #                      lower = 0,
  #                      weights = rep(1,3),
  #                      control = list(jemax = 2000))
  #   phi_hat <- opt3$coefficients
  # }
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2, unit should be the same as Aw
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2]/phi_hat[1] # Q_hat
  param[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
  param[4] <- phi_hat[3] * tau_hat # Aw_hat
  names(param) <- c("f0_hat", "Q_hat", "k_hat", "Aw_hat")
  return(param)
}
