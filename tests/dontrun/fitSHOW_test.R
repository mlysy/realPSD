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
                    bin_size = 100, method = c("lp", "mle", "nls"),
                    unit_conversion = FALSE) {
  # ---------- setup -----------
  method <- match.arg(method)
  Kb <- 1.381e-23           # Boltzmann's constant
  if(unit_conversion) {
    sig2 <- Kb*Temp/(k*pi*f0*Q) * 1e30 # variance, unit: fm2/Hz 
  } else {
    sig2 <- Kb*Temp/(k*pi*f0*Q) # variance, unit: m2/Hz
  }
  Rw <- Aw/sig2 # re-parameterization, note Rw is unitless 
  # phi <- c(f0, f0*Q, Rw) # parameter vector for SHOW model
  phi <- c(f0 + rnorm(1, 0, sqrt(f0)/10), 
    f0*Q + rnorm(1, 0, sqrt(f0*Q)/10), 
    Rw + rnorm(1,0,Rw/10)) 
  # psd values at each frequency point of f with given Q
  psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion) + Aw
  # generate the periodogram values
  Y <- sim_exp * psd * fs
  # convert Y to standard unit (otherwise the NLS optim would fail)
  if(unit_conversion) Y <- Y/1e30
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
  # opt <- optim(phi, fn = obj$fn, gr = obj$gr,
  #           method = "BFGS",
  #           control = list(maxit = 2000))
  if(method == "mle" || method == "lp"){
    opt <- optim(phi, fn = obj$fn, gr = obj$gr,
              method = "BFGS",
              control = list(maxit = 2000))
    phi_hat <- opt$par
  } else {
    # optimize Q (gamma), fix f0 and Rw
    opt1 <- pracma::lsqnonlin(fun = fn_fixed,
      x0 = phi, 
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)]) 
    # optimize Q and f0, fix Rw
    opt2 <- pracma::lsqnonlin(fun = fn_fixed, 
      x0 = opt1$x,
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = opt1$x[3])
    # optimize all three parameters
    opt3 <- pracma::lsqnonlin(fun = obj$fn, 
      x0 = opt2$x)
    # return phi_hat
    phi_hat <- opt3$x
  }
  # check convergence 
  if(opt$convergence != 0) 
    warning(paste0(method, " didn't converge!"))

  # output
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2, unit should be the same as Aw
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2]/phi_hat[1] # Q_hat
  param[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
  if(unit_conversion) {
    param[4] <- phi_hat[3] * tau_hat * 1e30 # Aw_hat, same unit
  } else {
    param[4] <- phi_hat[3] * tau_hat
  }
  names(param) <- c("f0_hat", "Q_hat", "k_hat", "Aw_hat")
  # return(param)

  # for matlab comparison
  fit_data <- rep(NA, 5)
  if(unit_conversion)
    fit_data[1] <- tau_hat * 1e30 # As
  else
    fit_data[1] <- tau_hat
  fit_data[2] <- param[2] # Q
  fit_data[3] <- param[1] # f0
  fit_data[4] <- param[4] # Aw
  fit_data[5] <- param[3] # k
  names(fit_data) <- c("As", "Q", "f0", "Aw", "k")
  return(fit_data)
}

# wrapper function of TMB ----- method 1 -----
MakePSDFun <- function(method = c("lp", "mle", "nls"), map, DLL) {
  if (method == "lp") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_nlp",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          map = map,
                          silent = TRUE, DLL = DLL)
  } else if (method == "nls") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_nlp",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          map = map,
                          silent = TRUE, DLL = DLL)
  } else if (method == "mle") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_nlp",
                                      f = matrix(fseq),
                                      Y = matrix(Y),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          map = map,
                          silent = TRUE, DLL = DLL)
  } else {
    stop("method should be chosen from lp, nls and mle.")
  }
  return(obj)
}
# wrapper of get_tau only be used after the last stage of optim
get_tau <- function(method = c("lp", "mle", "nls"), phi, DLL) {
  if(method == "lp"){
    get_tau <- function(phi) {
      gz <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_zeta",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = DLL)
      zeta <- gz$fn(phi)
      # correct the bias
      exp(zeta - bias)
    }
  } else if(method == "nls") {
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_tau",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = DLL)
      gt$fn(phi)
    }
  } else if (method == "mle") {
     get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_tau",
                                      f = matrix(fseq),
                                      Y = matrix(Y),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = DLL)
      gt$fn(phi)
    }
  } else {
    stop("method should be chosen from lp, nls and mle.")
  }
}

# wrapper functions ----- method 2 ------
#' @param obj TMB obj
#' @param theta Parameter vector
#' @param fixed_flag Vector of TRUE/FALSE indicating which dimension of theta should be fixed
#' @param fixed_phi Vector of fixed values, length(fixed_phi) == length(which(fixed_id == TRUE))
fn_fixed <- function(theta, obj, fixed_flag, fixed_phi) {
  # set space without chaning the original theta
  theta_full <- rep(NA, length(theta))
  # fix part of theta
  # theta_full[which(fixed_flag == TRUE)] <- fixed_phi
  theta_full[fixed_flag == TRUE] <- fixed_phi
  # fill the remaining part with theta
  # theta_full[which(fixed_flag != TRUE)] <- theta
  theta_full[!fixed_flag] <- theta[!fixed_flag]
  # return
  obj$fn(theta_full)
}
gr_fixed <- function(theta, obj, fixed_flag, fixed_phi) {
  theta_full <- rep(NA, length(theta))
  theta_full[fixed_flag == TRUE] <- fixed_phi
  theta_full[!fixed_flag] <- theta[!fixed_flag]
  obj$gr(theta_full)
}

