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
  phi <- c(f0, Q, Rw) # parameter vector for SHOW model
  # phi <- c(f0 + rnorm(1, 0, sqrt(f0)/10),
  #   Q + rnorm(1, 0, sqrt(Q)),
  #   Rw + rnorm(1,0, Rw/10))
  # psd values at each frequency point of f with given Q
  # psd <- psdSHO(fseq, f0, Q, k, Temp, unit_conversion) + Aw
  psd <- exp(logSHOW(fseq, f0, Q, tau = sig2, Rw))
  # generate the periodogram values
  Y <- sim_exp * psd * fs
  # convert Y to standard unit (otherwise the NLS optim would fail)
  if(unit_conversion) Y <- Y/1e30
  param <- fitSHOW_TMB(fseq, Y, bin_size, method, phi, Temp, Kb)
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
  c(obj$simulate(phi)$res)
}
# wrapper of nls_res with fixed parameters
nls_res_fixed <- function(phi, obj, fixed_flag, fixed_phi) {
  phi_full <- rep(NA, length(phi))
  phi_full[fixed_flag == TRUE] <- fixed_phi
  phi_full[!fixed_flag] <- phi[!fixed_flag]
  nls_res(phi_full, obj)
}
