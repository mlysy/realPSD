#--- fitting helper functions --------------------------------------------------

#' Efficiently calculate restricted frequency basis.
#'
#' @param frng Frequency range to restrict to.
#' @param fs Sampling frequency.
#' @param N Number of frequencies in the full basis.
#'
#' @return The restricted frequency basis, i.e., those in `0:N/N * fs/2` which are also inside `frng`.
get_fseq <- function(frng, fs, N) {
  # frequency range restricted to full basis
  df <- fs/(2*N) # frequency discretization
  seq(from = ceiling(frng[1]/df)*df, to = floor(frng[2]/df)*df, by = df)
}

#' Get `theta` from SHOW parameters.
#'
#' The parameter vector `theta` is defined as
#'
#' \preformatted{
#' theta = (f0, Q, Rw = Sw/tau, tau)`,
#' }
#' where `tau = KbT/(pi k f0 Q)`.
#'
#' @template param-k
#' @template param-f0
#' @template param-Q
#' @template param-Sw
#' @template param-Temp
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
#' @template param-Temp
#'
#' @return Numeric vector with named elements `(k, f0, Q, Sw)`.
get_par <- function(theta, Temp) {
  Kb <- 1.381e-23             # Boltzmann's constant
  Sw <- theta[3] * theta[4]
  k <- Kb*Temp/(theta[4]*pi*theta[1]*theta[2])
  setNames(c(k, theta[1], theta[2], Sw),
           nm = c("k", "f0", "Q", "Sw"))
}

#' Calculate bin size correction factor for `LP` estimator.
#'
#' @param B Bin size (integer).
bin_const <- function(B) log(B) - digamma(B)

#' Calculate SHOW PSD.
#'
#' @param fseq Frequency vector (units: Hz).
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw White noise power (m^2/Hz).
#'
#' @return The SHOW PSD evaluated at \code{fseq} (m^2/Hz).
show_psd <- function(fseq, k, f0, Q, Sw, Temp) {
  theta <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
  ## Kb <- 1.381e-23             # Boltzmann's constant
  ## tau <- Kb*Temp/(k*pi*f0*Q)
  ## Rw <- Sw/tau
  f2 <- (fseq/theta["f0"])^2
  theta["tau"] * (theta["Rw"] + 1/((1 - f2)^2 + f2/theta["Q"]^2))
}

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

show_fit_lp <- function(fseq, Ypsd, fs, Temp, bin_size, phi0,
                        fit_type = c("direct", "incremental"), ...) {
  fit_type <- match.arg(fit_type)
  fbar <- binning(fseq, bin_size = bin_size)
  Zbar <- log(binning(Ypsd, bin_size = bin_size))
  constZ <- mean(Zbar) # normalize to avoid numerical overflow
  obj <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                    method = "LP_nlp",
                                    fbar = as.matrix(fbar),
                                    Zbar = as.matrix(Zbar - constZ),
                                    fs = fs/exp(bin_const(bin_size))),
                        parameters = list(phi = as.matrix(c(0,0,0))),
                        silent = TRUE, DLL = "realPSD_TMBExports")
  exitflag <- NULL
  phi <- phi0
  if(fit_type == "incremental") {
    # fit f0 conditioned on everything else
    # don't use Brent because fn is unreliable far from mode...
    fixed <- c(FALSE, TRUE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", # lower = min(fbar), upper = max(fbar),
                 ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
    # fit f0 and Q conditioned on everything else
    fixed <- c(FALSE, FALSE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
  }
  if(all(exitflag == 0)) {
    # fit all three parameters at once
    fit <- optim(par = phi,
                 fn = obj$fn, gr = obj$gr, method = "BFGS", ...)
    phi <- fit$par
    exitflag <- c(exitflag, fit$convergence)
  }
  # construct final estimator
  theta <- c(exp(phi), tau = exp(obj$simulate(phi)$zeta + constZ))
  list(par = get_par(theta, Temp = Temp),
       value = obj$fn(phi),
       exitflag = exitflag)
}

show_fit_mle <- function(fseq, Ypsd, fs, Temp, phi0,
                        fit_type = c("direct", "incremental"), ...) {
  fit_type <- match.arg(fit_type)
  constY <- mean(Ypsd) # normalize to avoid numerical overflow
  obj <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                    method = "MLE_nlp",
                                    f = as.matrix(fseq),
                                    Y = as.matrix(Ypsd/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(c(0,0,0))),
                        silent = TRUE, DLL = "realPSD_TMBExports")
  exitflag <- NULL
  phi <- phi0
  if(fit_type == "incremental") {
    # fit f0 conditioned on everything else
    # don't use Brent because fn is unreliable far from mode...
    fixed <- c(FALSE, TRUE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", # lower = min(fbar), upper = max(fbar),
                 ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
    # fit f0 and Q conditioned on everything else
    fixed <- c(FALSE, FALSE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
  }
  if(all(exitflag == 0)) {
    # fit all three parameters at once
    fit <- optim(par = phi,
                 fn = obj$fn, gr = obj$gr, method = "BFGS", ...)
    phi <- fit$par
    exitflag <- c(exitflag, fit$convergence)
  }
  # construct final estimator
  theta <- c(exp(phi), tau = obj$simulate(phi)$tau * constY)
  list(par = get_par(theta, Temp = Temp),
       value = obj$fn(phi),
       exitflag = exitflag)
}

show_fit_nls <- function(fseq, Ypsd, fs, Temp, bin_size, phi0,
                        fit_type = c("direct", "incremental"), ...) {
  fit_type <- match.arg(fit_type)
  fbar <- binning(fseq, bin_size = bin_size)
  Ybar <- binning(Ypsd, bin_size = bin_size)
  constY <- mean(Ybar) # normalize to avoid numerical overflow
  obj <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                    method = "NLS_nlp",
                                    fbar = as.matrix(fbar),
                                    Ybar = as.matrix(Ybar/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(c(0,0,0))),
                        silent = TRUE, DLL = "realPSD_TMBExports")
  exitflag <- NULL
  phi <- phi0
  if(fit_type == "incremental") {
    # fit f0 conditioned on everything else
    # don't use Brent because fn is unreliable far from mode...
    fixed <- c(FALSE, TRUE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", # lower = min(fbar), upper = max(fbar),
                 ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
    # fit f0 and Q conditioned on everything else
    fixed <- c(FALSE, FALSE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
  }
  if(all(exitflag == 0)) {
    # fit all three parameters at once
    fit <- optim(par = phi,
                 fn = obj$fn, gr = obj$gr, method = "BFGS", ...)
    phi <- fit$par
    exitflag <- c(exitflag, fit$convergence)
  }
  # construct final estimator
  theta <- c(exp(phi), tau = obj$simulate(phi)$tau * constY)
  list(par = get_par(theta, Temp = Temp),
       value = obj$fn(phi),
       exitflag = exitflag)
}

#' Displays error message and returns NA.
#'
#' To be used as `error` argument of  `tryCatch` function.
#' @param err Error message.
catch_error <- function(err) {
  print(err)
  NA
}
