#' SHOWF fit method non-linear least squares (NLS)
#' @param fseq Frequency vector (units: Hz).
#' @param Ypsd PSD vector.
#' @param fs Sampling frequency.
#' @param Temp True temperature.
#' @param bin_size Bin size.
#' @param bin_type Either "mean" or "median".  The former is more efficient, the latter is more robust, so better for preliminary SHOW fit to denoise against.
#' @param phi0 Parameter value (transformed scale) to initialize optimization.
#' @param fit_type If "direct", fit all parameters at once.  If "incremental", do one, then two, then three, etc.
#' @param optimizer Either "optim" (in R) or "Adam" (supplied by realPSD). For now, Adam's hyperparameter tunning is not fully supported. Learning rate and nsteps are preset. `Adam` is only used for `direct` fitting.
#' @param vcov If TRUE, return the numerical Hessian matrix of the original SHOF parameter k, f0, Q, Af
#' @param ... Additional arguments to [stats::optim()].
#' 
#' @return A list with elements 
#' \describe{
#'   \item{`par`}{The fitted parameter value which is the minimizer of the objective function.}
#'   \item{`Temp`}{The temperature constant used in calculation, Kelvin.}
#'   \item{`value`}{The value of the objective function (negative loglikelihood) at the given `par`.}
#'   \item{`cov`}{The numerical variance-covariance matrix for all the parameters (original scale).}
#'   \item{`jac`}{The numerical Jacobian matrix of the SHOF model parameter k, f0, Q, Af.}
#'   \item{`exitflag`}{The exit flag given by the optimizer, for debug purposes}
#' }
showf_fit_nls <- function(fseq, Ypsd, fs, Temp,
                         bin_size, bin_type, phi0,
                        fit_type = c("direct", "incremental"),
                        optimizer = c("optim", "Adam"),
                        vcov = FALSE, get_jac = FALSE, ...) {
  fit_type <- match.arg(fit_type)
  optimizer <- match.arg(optimizer)
  fbar <- binning(fseq, bin_size = bin_size, bin_type = bin_type)
  Ybar <- binning(Ypsd, bin_size = bin_size, bin_type = bin_type)
  constY <- mean(Ybar) # normalize to avoid numerical overflow
  map <- list(as.factor(c(1,2,NA,4,5)))
  obj <- TMB::MakeADFun(data = list(model = "SHOWF_log",
                                    method = "NLS_nlp",
                                    fbar = as.matrix(fbar),
                                    Ybar = as.matrix(Ybar/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(rep(0,5))),
                        map = map,
                        silent = TRUE, DLL = "realPSD_TMBExports")
  exitflag <- NULL
  phi <- phi0
  if(fit_type == "incremental") {
    # fit f0 conditioned on everything else
    # don't use Brent because fn is unreliable far from mode...
    fixed <- c(FALSE, TRUE, TRUE, TRUE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", # lower = min(fbar), upper = max(fbar),
                 ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
    # fit f0 and Q conditioned on everything else
    fixed <- c(FALSE, FALSE, TRUE, TRUE, TRUE)
    fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", ...)
    phi[!fixed] <- fit$par
    exitflag <- c(exitflag, fit$convergence)
    if(any(exitflag) != 0) break
    # fit f0, Q, Rf conditioned on alpha
    fixed <- c(FALSE, FALSE, TRUE, FALSE, TRUE)
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
    if(optimizer == "optim") {
      fixed <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
      # fit <- optim(par = phi,
      #            fn = obj$fn, gr = obj$gr, method = "BFGS", ...)
      fit <- optim(par = phi[!fixed],
                 fn = fn_fixed, gr = gr_fixed,
                 obj = obj, fixed = fixed, phi0 = phi,
                 method = "BFGS", ...)
      phi[!fixed] <- fit$par
    } else if(optimizer == "Adam") {
      fit <- adam(theta0 = phi, fn = obj$fn, gr = obj$gr, nsteps = 300,
                alpha = 1e-4, ...)
      phi <- fit$par
    }
    exitflag <- c(exitflag, fit$convergence)
  }
  # construct final estimator
  theta <- c(exp(phi), tau = obj$simulate(phi)$tau * constY)
  theta <- setNames(theta, nm = c("f0", "Q", "Rw", "Rf", "alpha", "tau"))
  par_opt <- get_par_showf(theta, Temp = Temp)
  # hessian
  he <- NULL
  if(vcov) {    
    obj_nll <- TMB::MakeADFun(data = list(model = "SHOWF_log",
                                    method = "NLS_nll",
                                    fbar = as.matrix(fseq),
                                    Ybar = as.matrix(Ypsd/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(rep(0,5)), tau = 0),
                        map = map,
                        silent = TRUE, DLL = "realPSD_TMBExports")
    he <- numDeriv::hessian(func = function(par) {
      # convert original scale (par) to computational scale (phi & zeta)
      phi_tau <- get_phi(par, Temp = Temp, method = "NLS", model = "SHOWF", const = constY)
      # feed these into the negative loglikelihood on the computational scale
      obj_nll$fn(phi_tau)
    }, x = par_opt) 
    he <- he[c(1:3,6), c(1:3,6)]
    cov <- chol2inv(chol(he)) 
  }
  # Jacobian
  jac <- NULL
  if(get_jac) {
    obj_res <- TMB::MakeADFun(data = list(model = "SHOWF_log",
                                    method = "NLS_res",
                                    fbar = as.matrix(fbar),
                                    Ybar = as.matrix(Ybar/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(rep(0,5))),
                        map = map,
                        silent = TRUE,
                        ADreport = TRUE,
                        DLL = "realPSD_TMBExports")
    nls_res2 <- function(par) {
      phi_tau <- get_phi(par, Temp = Temp, method = "NLS", model = "SHOWF", const = constY)
      (obj_res$fn(phi_tau[1:5]))^2
    }
    jac <- numDeriv::jacobian(nls_res2, x = par_opt)
  }
  list(par = get_par_showf(theta, Temp = Temp),
       value = obj$fn(phi),
       cov = cov, jac = jac,
       exitflag = exitflag)
}
