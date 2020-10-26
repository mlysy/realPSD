# ---------- internal helper functions used by show fit methods show_fit_*.R ----------
#' ADAM optimizer.
#'
#' @param theta0 Initial parameter value (numeric vector).
#' @param fn Objective function to minimize.
#' @param gr Gradient function of `fn`.
#' @param nsteps Number of optimization steps to perform.
#' @param alpha,beta,eps Algorithm tuning parameters.  `alpha` and `eps` are scalars, `beta` is a numeric vector of length 2.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`par`}{The parameter value at the minimum of the `nsteps` optimization steps.}
#'   \item{`objective`}{The minimum value of the objective function.}
#'   \item{`Theta`}{A matrix where each row is the parameter value at the given step.}
#'   \item{`Fn`}{A vector where each element is the objective function value at the given step.}
#' }
adam <- function(theta0, fn, gr, nsteps,
                 alpha = 1e-3, beta = c(.9, .999), eps = 1e-8) {
  Fn <- rep(NA, nsteps)
  Theta <- matrix(NA, nsteps, length(theta0))
  colnames(Theta) <- names(theta0)
  mt <- 0
  vt <- 0
  thetat <- theta0
  betat <- c(1,1)
  for(tt in 1:nsteps) {
    gt <- gr(thetat)
    mt <- beta[1] * mt + (1-beta[1]) * gt
    vt <- beta[2] * vt + (1-beta[2]) * gt^2
    betat <- betat * beta
    mhat <- mt/(1-betat[1])
    vhat <- vt/(1-betat[2])
    thetat <- thetat - alpha * mhat / (sqrt(vhat) + eps)
    Theta[tt,] <- thetat
    Fn[tt] <- fn(thetat)
  }
  imin <- which.min(Fn)
  list(par = Theta[imin,], objective = Fn[imin],
       Theta = Theta, Fn = Fn)
}

#' Evaluate a \pkg{TMB} objective function at fixed parameter values.
#'
fn_fixed <- function(phi, obj, fixed, phi0) {
  x <- phi0
  x[!fixed] <- phi
  obj$fn(x)
}

#' Evaluate a \pkg{TMB} gradient function at fixed parameter values.
#'
gr_fixed <- function(phi, obj, fixed, phi0) {
  x <- phi0
  x[!fixed] <- phi
  obj$gr(x)[!fixed]
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

#' Recover the SHOWF parameters from `theta`.
#'
#' @param theta Parameter vector (f0, Q, Rw, Rf, alpha, tau).
#' @param Temp Temperature (Kelvin).
#'
#' @return Numeric vector with named elements `(k, f0, Q, Sw, Af, alpha)`.
get_par_showf <- function(theta, Temp) {
  Kb <- 1.381e-23             # Boltzmann's constant
  Sw <- theta[3] * theta[6] # Rw * tau
  k <- Kb*Temp/(theta[6]*pi*theta[1]*theta[2])
  Af <- theta[4] * theta[6] # Rf * tau
  alpha <- theta[5]
  setNames(c(k, theta[1], theta[2], Sw, Af, alpha),
           nm = c("k", "f0", "Q", "Sw", "Af", "alpha"))
}

#' Recover the SHOF parameters from `theta`.
#'
#' @param theta Parameter vector (f0, Q, Rf, alpha, tau).
#' @param Temp Temperature (Kelvin).
#'
#' @return Numeric vector with named elements `(k, f0, Q, Af, alpha)`.
get_par_shof <- function(theta, Temp) {
  Kb <- 1.381e-23             # Boltzmann's constant
  k <- Kb*Temp/(theta[5]*pi*theta[1]*theta[2])
  Af <- theta[3] * theta[5] # Rf * tau
  alpha <- theta[4]
  setNames(c(k, theta[1], theta[2], Af, alpha),
           nm = c("k", "f0", "Q", "Af", "alpha"))
}

# #' Calculate bin size correction factor for `LP` estimator.
# #'
# #' @param B Bin size (integer).
# #'
# #' @return The bin size correction factor `log(B) - digamma(B)`.
# #' @note This doesn't give the correct value for median estimator...
# bin_const <- function(B) log(B) - digamma(B)
