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
