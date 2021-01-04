#' Convert SHOW parameterization to natural parameterization.
#' @param par SHOW parameters `(k, f0, Q, Sw)` or `(k, f0, Q, Sw, Af, alpha)` or `(k, f0, Q, Af, alpha)`.
#' @param method Fitting method
#' @param model either "SHOW" or "SHOWF", if model == "SHOW", par = c(k,f0,Q,Sw); if model == "SHOWF", par = c(k,f0,Q,Sw,Af,alpha)
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
#' @export
get_phi <- function(par, method = c("MLE", "NLS", "LP"), model = c("SHOW", "SHOWF", "SHOF"), Temp, const) {
  Kb <- 1.381e-23
  k <- par[1]
  f0 <- par[2]
  Q <- par[3]
  Sw <- par[4]
  tau <- Kb*Temp/(k*pi*f0*Q)
  Rw <- Sw/tau
  if(model == "SHOW") {
    # normalized tau based on fitting method to avoid numerical overflow
    if(method == "MLE") {
      tau <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), tau), nm=c("log f0", "log Q", "log Rw", "tau"))
      # phi_zeta <- setNames(c(f0, Q, Rw, tau), nm=c("f0", "Q", "Rw", "tau"))
    } else if (method == "NLS") {
      tau <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), tau), nm=c("log f0", "log Q", "log Rw", "tau"))
      # phi_zeta <- setNames(c(f0, Q, Rw, tau), nm=c("f0", "Q", "Rw", "tau"))
    } else if (method == "LP") {
      zeta <- log(tau) - const
      # zeta <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), zeta), nm=c("log f0", "log Q", "log Rw", "zeta"))
      # phi_zeta <- setNames(c(f0, Q, Rw, zeta), nm=c("f0", "Q", "Rw", "zeta"))
    }
  } else if (model == "SHOWF") {
    Af <- par[5]
    alpha <- par[6]
    Rf <- Af/tau
    # normalized tau based on fitting method to avoid numerical overflow
    if(method == "MLE") {
      tau <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), log(Rf), log(alpha), tau), nm=c("log f0", "log Q", "log Rw", "log Rf", "log alpha", "tau"))
    } else if (method == "NLS") {
      tau <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), log(Rf), log(alpha), tau), nm=c("log f0", "log Q", "log Rw", "log Rf", "log alpha", "tau"))
    } else if (method == "LP") {
      zeta <- log(tau) - const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rw), log(Rf), log(alpha), zeta), nm=c("log f0", "log Q", "log Rw", "log Rf", "log alpha", "zeta"))
    }
  } else if (model == "SHOF") {
    # par = c(k,f0,Q,Af,alpha)
    Af <- par[4]
    alpha <- par[5]
    Rf <- Af/tau
    if(method == "MLE") {
      tau <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rf), log(alpha), tau), nm=c("log f0", "log Q", "log Rf", "log alpha", "tau"))
    } else if (method == "NLS") {
      tau <- tau/const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rf), log(alpha), tau), nm=c("log f0", "log Q", "log Rf", "log alpha", "tau"))
    } else if (method == "LP") {
      zeta <- log(tau) - const
      phi_zeta <- setNames(c(log(f0), log(Q), log(Rf), log(alpha), zeta), nm=c("log f0", "log Q", "log Rf", "log alpha", "zeta"))
    }
  }
  return(phi_zeta)
}
