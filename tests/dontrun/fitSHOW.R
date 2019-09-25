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
#' @param unit_conversion Logical, if TRUE, use fm2/Hz instead of m2/Hz
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
  Rw <- Aw/sig2 # re-parameterization, note we input Aw with unit fm2/Hz
  phi <- c(f0, f0*Q, Rw) # parameter vector for SHOW model
  # phi <- c(f0 + rnorm(1, 0, sqrt(f0)/10), 
  #   f0*Q + rnorm(1, 0, sqrt(f0*Q)/10), 
  #   Rw + rnorm(1,0, Rw/10)) 
  # psd values at each frequency point of f with given Q
  psd <- psdSHO(fseq, f0, Q, k, Temp, unit_conversion) + Aw
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
  # phi_hat <- opt$par
  if(method == "mle" || method == "lp"){
    # ---------- optim all-at-once ---------
    # opt <- optim(phi, fn = obj$fn, gr = obj$gr,
    #           method = "BFGS",
    #           control = list(maxit = 2000))
    # # check convergence for MLE and LP
    # if(opt$convergence != 0) 
    #   warning(paste0(method, " didn't converge!"))
    # phi_hat <- opt$par
    # ----- optim step-by-step -----
    start <- phi
    opt1 <- optim(start, fn = fn_fixed, 
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)],
      method = "Nelder-Mead") 
    start[2] <- opt1$par[2]
    opt2 <- optim(start, fn = fn_fixed,
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3],
      method = "Nelder-Mead")
    start[c(1,2)] <- opt2$par[c(1,2)]
    # opt3 <- optim(start, fn = obj$fn)
    opt3 <- optim(start, fn = obj$fn, gr = obj$gr, method = "BFGS")
    phi_hat <- opt3$par
    # ----- fminsearch ------
    # start <- phi
    # opt1 <- pracma::fminsearch(fn = fn_fixed,
    #   x0 = start,
    #   obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)],
    #   method = "Nelder-Mead") 
    # start[2] <- opt1$xmin[2]
    # opt2 <- pracma::fminsearch(fn = fn_fixed,
    #   x0 = start,
    #   obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3],
    #   method = "Nelder-Mead")
    # start[c(1,2)] <- opt2$xmin[c(1,2)]
    # opt3 <- pracma::fminsearch(fn = obj$fn, 
    #   x0 = start, method = "Nelder-Mead")
    # phi_hat <- opt3$xmin
  } else {
    # ---------- lsqnonlin ---------
    # set some option parameters to avoid errors
    # tolx <- 1/(sum(phi^2)) # in order to compensate the squared norm of the parameter vector
    # tolg <- max(10*abs(obj$gr(phi)))
    tolx <- .Machine$double.eps
    tolg <- .Machine$double.eps
    # optimize Q (gamma), fix f0 and Rw
    opt1 <- pracma::lsqnonlin(fun = nls_res_fixed,
      x0 = phi, 
      options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)]) 
    # if(opt1$errno != 1) warning("NLS didn't converge at step 1.")
    # optimize Q and f0, fix Rw
    opt2 <- pracma::lsqnonlin(fun = nls_res_fixed, 
      x0 = opt1$x,
      options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3])
    # if(opt1$errno != 1) warning("NLS didn't converge at step 2.")
    # optimize all three parameters
    opt3 <- pracma::lsqnonlin(fun = nls_res_fixed, 
      x0 = opt2$x,
      options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
      obj = obj, fixed_flag = c(0,0,0), fixed_phi = NULL)
    # if(opt1$errno != 1) warning("NLS didn't converge at step 3.")
    # return phi_hat
    phi_hat <- opt3$x
    # ---------- minpack.lm::nls.lm ---------
    # opt1 <- minpack.lm::nls.lm(par = phi, lower = rep(0,3), 
    #   fn = nls_res_fixed,
    #   obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)]) 
    # opt2 <- minpack.lm::nls.lm(par = opt1$par, lower = rep(0,3), 
    #   fn = nls_res_fixed,
    #   obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3])  
    # opt3 <- minpack.lm::nls.lm(par = opt2$par, lower = rep(0,3), 
    #   fn = nls_res_fixed,
    #   obj = obj, fixed_flag = c(0,0,0), fixed_phi = NULL)  
    # phi_hat <- opt3$par
  }
  # output
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2, unit should be the same as Aw
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2]/phi_hat[1] # Q_hat
  param[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
  if(unit_conversion) {
    # param[4] <- phi_hat[3] * tau_hat * 1e30 # Aw_hat, same unit
    param[4] <- Kb * Temp / (pi * phi_hat[2]) * 1e30
  } else {
    # param[4] <- phi_hat[3] * tau_hat
    param[4] <- Kb * Temp / (pi * phi_hat[2])
  }
  names(param) <- c("f0_hat", "Q_hat", "k_hat", "Aw_hat")
  return(param)
}

# wrapper functions ----- method 1 ------
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
# wrapper function of the vector of residuals for NLS
nls_res <- function(phi, obj) {
  c(obj$simulate(phi)$RES)
} 
# wrapper of nls_res with fixed parameters
nls_res_fixed <- function(phi, obj, fixed_flag, fixed_phi) {
  phi_full <- rep(NA, length(phi))
  phi_full[fixed_flag == TRUE] <- fixed_phi
  phi_full[!fixed_flag] <- phi[!fixed_flag]
  nls_res(phi_full, obj)
}
