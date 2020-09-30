#' Fit SHOW model
#' @param fseq Frequency vector (units: Hz).
#' @param Ypsd PSD vector, returned by realPSD::periodogram() which is |fft(X_t)|^2/(N*fs).
#' @param fs Sampling frequency.
#' @param Temp True temperature.
#' @param method One of the three fitting methods "LP", "NLS" and "MLE".
#' @param bin_size Bin size.
#' @param bin_type Either "mean" or "median".  The former is more efficient, the latter is more robust, so better for preliminary SHOW fit to denoise against.
#' @param phi0 Parameter value (transformed scale) to initialize optimization.
#' @param fit_type If "direct", fit all parameters at once.  If "incremental", do one, then two, then three, etc.
#' @param getHessian If TRUE, return the numerical Hessian matrix of the original SHOW parameters f0, Q, Rw
#' @param getJacobian If TRUE, return the numerical Jacobian of the SHOW parameters f0, Q, Rw, only valid if method == "NLS".
#' @param ... Additional arguments to [stats::optim()].
#'
#' @return A list with elements 
#' \describe{
#'   \item{`par`}{The fitted parameter value which is the minimizer of the objective function.}
#'   \item{`Temp`}{The temperature constant used in calculation, Kelvin.}
#'   \item{`value`}{The value of the objective function (negative loglikelihood) at the given `par`.}
#'   \item{`he`}{The numerical Hessian matrix of the SHOW model parameter f0, Q, Rw.}
#'   \item{`jac`}{The numerical Jacobian matrix of the SHOW model parameter f0, Q, Rw, only vailable if method == "NLS".}
#'   \item{`exitflag`}{The exit flag given by the optimizer, for debug purposes}
#' }
#' 
#' @note Ypsd used as input argument should be returned by realPSD::periodogram which is 1/(N*fs) * |fft(X_t)|^2, but TMB method to fit the SHOW model assumes unscaled periodogram, so we first transform Ypsd in the function. The TMB fitting methods should be modified later to account for this.
#' 
#' @export
show_fit <- function(fseq, Ypsd, fs, Temp,
                     method = c("lp", "nls", "mle"),
                     bin_size, bin_type, phi0,
                     fit_type = c("direct", "incremental"),
                     getHessian = FALSE,
                     getJacobian = FALSE, ...) {
  Ypsd <- Ypsd * fs # since Ypsd is returned by periodogram which has been scaled by 1/fs.
  method <- match.arg(method)
  if(method == "lp") {
    out <- show_fit_lp(fseq = fseq, Ypsd = Ypsd, fs = fs, Temp = Temp,
                       bin_size = bin_size, bin_type = bin_type,
                       phi0 = phi0, fit_type = fit_type, 
                       getHessian = getHessian, ...)
  } else if(method == "nls") {
    out <- show_fit_nls(fseq = fseq, Ypsd = Ypsd, fs = fs, Temp = Temp,
                        bin_size = bin_size, bin_type = bin_type,
                        phi0 = phi0, fit_type = fit_type, 
                        getHessian = getHessian, ...)
  } else if(method == "mle") {
    out <- show_fit_mle(fseq = fseq, Ypsd = Ypsd, fs = fs, Temp = Temp,
                        phi0 = phi0, fit_type = fit_type, 
                        getHessian = getHessian, getJacobian = getJacobian, ...)
  }
  out
}

#--- fitting helper functions --------------------------------------------------
#' Get `theta` from SHOW parameters.
#'
#' The parameter vector `theta` is defined as
#' ```
#' theta = (f0, Q, Rw = Sw/tau, tau),
#' ```
#' where `tau = KbT/(pi k f0 Q)`.
#'
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw White noise power (m^2/Hz).
#' @param Temp Temperature (Kelvin).
#'
#' @return Numeric vector with named elements `(f0, Q, Rw, tau)`.
get_theta <- function(k, f0, Q, Sw, Temp) {
  Kb <- 1.381e-23             # Boltzmann's constant
  tau <- Kb*Temp/(k*pi*f0*Q)
  Rw <- Sw/tau
  setNames(c(f0, Q, Rw, tau), nm = c("f0", "Q", "Rw", "tau"))
}

#' Recover the SHOW parameters from `theta`.
#'
#' @param theta Parameter vector as defined by `get_theta`.
#' @param Temp Temperature (Kelvin).
#'
#' @return Numeric vector with named elements `(k, f0, Q, Sw)`.
get_par <- function(theta, Temp) {
  Kb <- 1.381e-23             # Boltzmann's constant
  Sw <- theta[3] * theta[4]
  k <- Kb*Temp/(theta[4]*pi*theta[1]*theta[2])
  setNames(c(k, theta[1], theta[2], Sw),
           nm = c("k", "f0", "Q", "Sw"))
}

# #' Calculate bin size correction factor for `LP` estimator.
# #'
# #' @param B Bin size (integer).
# #'
# #' @return The bin size correction factor `log(B) - digamma(B)`.
# #' @note This doesn't give the correct value for median estimator...
# bin_const <- function(B) log(B) - digamma(B)

#' Calculate SHOW PSD.
#'
#' @param fseq Frequency vector (units: Hz).
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw White noise power (m^2/Hz).
#'
#' @return The SHOW PSD evaluated at `fseq` (m^2/Hz).
show_psd <- function(fseq, k, f0, Q, Sw, Temp) {
  theta <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
  ## Kb <- 1.381e-23             # Boltzmann's constant
  ## tau <- Kb*Temp/(k*pi*f0*Q)
  ## Rw <- Sw/tau
  f2 <- (fseq/theta["f0"])^2 
  theta["tau"] * (theta["Rw"] + 1/((1 - f2)^2 + f2/theta["Q"]^2))
}


#' Convert SHOW parameterization to natural parameterization.
#' @param par SHOW parameters `(k, f0, Q, Sw)`.
#' @param method Fitting method
#' @param Temp Temperature constant
#' @param const Normalizing constant applied to both PSD inputs and `tau` in each fitting method (see 'Notes').
#' @return Numeric vector with elements that you can plug back into TMB methods `LP_nll`, `NLS_nll`, or `MLE_nll`.
#'
#' @details Since we want the result in the original basis, let's compute the hessian numerically.  That is, if `par_hat` is the mode and let's say we're doing LP, the hessian is obtained with
#' ```
#' # assume we are inside show_fit_*,
#' # such that `bin_size`, `constZ`, `Zbar`, etc. have been defined
#' # (and properly normalized).
#' numDeriv::hessian(
#'             x = par_hat,
#'             func = function(par) {
#'               phi_zeta <- get_phi(par, method = "LP", const = constZ)
#'               obj <- TMB::MakeADFun(data = list(
#'                                       model = "SHOW_log",
#'                                       method = "LP_nll",
#'                                       fbar = as.matrix(fbar),
#'                                       Zbar = as.matrix(Zbar - constZ),
#'                                       fs = fs/exp(bin_const(bin_size))),
#'                                     parameters = list(phi = as.matrix(c(0,0,0)),
#'                                                       zeta = 0),
#'                                     silent = TRUE, DLL = "realPSD_TMBExports")
#'               obj$fn(phi_zeta)
#'             })
#'
#' ```
#'
#' @note
#' - [show_fit_mle()] and [show_fit_nls()] normalize the inputs to avoid numerical overflow.  This should be done here as well, which means there's a conversion factor to apply to `tau`.  My suggestion is to add the standard error calculation directly to the `show_fit_*` functions.
#' - `show_fit_*` have something called `constZ` or `constY` which is used to rescale the PSD inputs to avoid numerical overflow.  The scale factor `tau` gets corrected accordingly, and the same thing should be done here as well with the parameter `const`.
get_phi <- function(par, method = c("MLE", "NLS", "LP"), Temp, const) {
  Kb <- 1.381e-23
  k <- par[1]
  f0 <- par[2]
  Q <- par[3]
  Sw <- par[4]
  tau <- Kb*Temp/(k*pi*f0*Q)
  Rw <- Sw/tau
  # normalized tau based on fitting method to avoid numerical overflow
  if(method == "MLE") {
    tau <- tau/const
    phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), tau), nm=c("f0", "Q", "Rw", "tau"))
    # phi_zeta <- setNames(c(f0, Q, Rw, tau), nm=c("f0", "Q", "Rw", "tau"))
  } else if (method == "NLS") {
    tau <- tau/const
    phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), tau), nm=c("f0", "Q", "Rw", "tau"))
    # phi_zeta <- setNames(c(f0, Q, Rw, tau), nm=c("f0", "Q", "Rw", "tau"))
  } else if (method == "LP") {
    zeta <- log(tau) - const
    # zeta <- tau/const
    phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), zeta), nm=c("f0", "Q", "Rw", "zeta"))
    # phi_zeta <- setNames(c(f0, Q, Rw, zeta), nm=c("f0", "Q", "Rw", "zeta"))
  }
  return(phi_zeta)
}
