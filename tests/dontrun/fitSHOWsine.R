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
#' @param bin_size Integer number, bin size.
#' @param method Fitting method, i.e. lp, mle, nls.
#' @param unit_conversion Logical, if TRUE, use fm2/Hz instead of m2/Hz
fitSHOWsine <- function(fseq, sim_cnorm, f0, fs, Q, k, Temp, Aw,
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
  phi <- c(f0, f0*Q, Rw) # parameter vector for SHOW model
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
  sin_fft <- fft_sin(fseq, f0, Q, fs, unit_conversion)
  Y <- sim_cnorm * sqrt(psd * fs)
  Y <- (Y + sin_fft) * Conj(Y + sin_fft)
  Y <- Re(Y)
  # remove sine wave noise 
  if(remove_noise) {
    freq_range <- c(f0-f0/sqrt(2), f0+f0/sqrt(2))
    Y <- psd_denoise(fseq, psd_noise = Y, 
      Q, f0, k, Temp, unit_conversion, Aw, freq_range)
  }
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
  if(method == "mle" || method == "lp"){
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
    opt3 <- optim(start, fn = obj$fn, gr = obj$gr, 
      method = "L-BFGS-B", lower = rep(0,3))
    phi_hat <- opt3$par
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
    # opt3 <- pracma::lsqnonlin(fun = nls_res_fixed, 
    #   x0 = opt2$x,
    #   options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
    #   obj = obj, fixed_flag = c(0,0,0), fixed_phi = NULL)
    opt3 <- optim(par = opt2$x, fn = obj$fn, gr = obj$gr, 
      method = "L-BFGS-B", lower = rep(0,3))
    # if(opt1$errno != 1) warning("NLS didn't converge at step 3.")
    # return phi_hat
    phi_hat <- opt3$par
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
# # FFT (discrete Fourier transform) of sine wave noise
# fft_sin <- function(fseq, f0, Q, fs, unit_conversion) {
#   N <- length(fseq)
#   Const <- 1e30                # unit conversion, 1 m2 = 1e30 fm2
#   if(!unit_conversion) {
#     D <- (Q^0.5) * 3.5e3 / Const
#   } else {
#     D <- (Q^0.5) * 3.5e3 # 3.5e3 is copied from the MATLAB code
#   }
#   dT <- 1/fs
#   xi <- rnorm(1, f0, 10)
#   phi <- runif(1, 0, 2*pi)
#   sin_fft <- D/(2*1i * sqrt(N)) * (
#     exp(phi*1i) * (exp(2*pi*1i*(xi-fseq)*dT*N)-1)/(exp(2*pi*1i*(xi-fseq)*dT)-1) -
#     exp(-phi*1i) * (exp(-2*pi*1i*(xi+fseq)*dT*N)-1)/(exp(-2*pi*1i*(xi+fseq)*dT)-1)
#   )
#   return(sin_fft)
# }
