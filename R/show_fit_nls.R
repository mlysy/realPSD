#' SHOW fit method non-linear least squares (NLS)
#' @param fseq Frequency vector (units: Hz).
#' @param Ypsd PSD vector.
#' @param fs Sampling frequency.
#' @param Temp True temperature.
#' @param bin_size Bin size.
#' @param bin_type Either "mean" or "median".  The former is more efficient, the latter is more robust, so better for preliminary SHOW fit to denoise against.
#' @param phi0 Parameter value (transformed scale) to initialize optimization.
#' @param fit_type If "direct", fit all parameters at once.  If "incremental", do one, then two, then three, etc.
#' @param optimizer Either "optim" (in R) or "Adam" (supplied by realPSD). For now, Adam's hyperparameter tunning is not fully supported. Learning rate and nsteps are preset. `Adam` is only used for `direct` fitting.
#' @param vcov If TRUE, return the numerical variance-covariance matrix for all the parameters (original scale)
#' @param get_jac If TRUE, return the numerical Jacobian of the SHOW parameters f0, Q, Rw, only valid if method == "NLS".
#' @param ... Additional arguments to [stats::optim()].
#' 
#' @return A list with elements 
#' \describe{
#'   \item{`par`}{The fitted parameter value which is the minimizer of the objective function.}
#'   \item{`Temp`}{The temperature constant used in calculation, Kelvin.}
#'   \item{`value`}{The value of the objective function (negative loglikelihood) at the given `par`.}
#'   \item{`cov`}{The numerical variance-covariance matrix for all the parameters (original scale).}
#'   \item{`jac`}{The numerical Jacobian matrix of the SHOW model parameter f0, Q, Rw.}
#'   \item{`exitflag`}{The exit flag given by the optimizer, for debug purposes}
#' }
show_fit_nls <- function(fseq, Ypsd, fs, Temp,
                         bin_size, bin_type, phi0,
                        fit_type = c("direct", "incremental"),
                        optimizer = c("optim", "Adam"),
                        vcov = FALSE, get_jac = FALSE, ...) {
  fit_type <- match.arg(fit_type)
  optimizer <- match.arg(optimizer)
  fbar <- binning(fseq, bin_size = bin_size, bin_type = bin_type)
  Ybar <- binning(Ypsd, bin_size = bin_size, bin_type = bin_type)
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
    if(optimizer == "optim") {
      fit <- optim(par = phi,
                 fn = obj$fn, gr = obj$gr, method = "BFGS", ...)
    } else if(optimizer == "Adam") {
      fit <- adam(theta0 = phi, fn = obj$fn, gr = obj$gr, nsteps = 300,
                alpha = 1e-4, ...)
    }
    phi <- fit$par
    exitflag <- c(exitflag, fit$convergence)
  }
  # construct final estimator
  theta <- c(exp(phi), tau = obj$simulate(phi)$tau * constY)
  theta <- setNames(theta, nm = c("f0", "Q", "Rw", "tau"))
  par_opt <- get_par(theta, Temp = Temp)
  # hessian
  cov <- NULL
  if(vcov) {    
    obj_nll <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                    method = "NLS_nll",
                                    fbar = as.matrix(fseq),
                                    Ybar = as.matrix(Ybar/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(c(0,0,0)), tau = 0),
                        silent = TRUE, DLL = "realPSD_TMBExports")
    # phi_tau <- c(exp(phi), obj$simulate(phi)$tau)
    # he <- numDeriv::hessian(func = obj_nll$fn, x = phi_tau)
    # he <- he[1:3, 1:3] # truncate the row and col wrt tau
    he <- numDeriv::hessian(func = function(par) {
      # convert original scale (par) to computational scale (phi & zeta)
      phi_tau <- get_phi(par, Temp = Temp, method = "NLS", model = "SHOW", const = constY)
      # feed these into the negative loglikelihood on the computational scale
      obj_nll$fn(phi_tau)
    }, x = par_opt, method.args = list(zero.tol = .Machine$double.eps, r=6)) # we need to set a smaller zero.tol otherwise NaN will be produced
    cov <- chol2inv(chol(he)) 
  }
  # Jacobian
  jac <- NULL
  if(get_jac) {
    # obj_res <- TMB::MakeADFun(data = list(model = "SHOW_log",
    #                                 method = "NLS_res",
    #                                 fbar = as.matrix(fbar),
    #                                 Ybar = as.matrix(Ybar/constY),
    #                                 fs = fs),
    #                     parameters = list(phi = as.matrix(c(0,0,0))),
    #                     silent = TRUE,
    #                     ADreport = TRUE,
    #                     DLL = "realPSD_TMBExports")
    # nls_res2 <- function(par) {
    #   phi_tau <- get_phi(par, Temp = Temp, method = "NLS", model = "SHOW", const = constY)
    #   setNames((obj_res$fn(phi_tau[1:3]))^2, nm = NULL)
    # }
    obj_ufun <- TMB::MakeADFun(data = list(model = "SHOW_log", 
                                      method = "UFun", 
                                      f = matrix(fbar),
                                      fs = fs),
                           parameters = list(phi = matrix(phi0)),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nls_res2 <- function(par) {
      phi_tau <- get_phi(par, Temp = Temp, method = "NLS", model = "SHOW", const = constY)
      ufun_vec <- fs * obj_ufun$simulate(phi)$U
      (Ybar/constY - phi_tau[5] * ufun_vec)^2
    }
    jac <- numDeriv::jacobian(nls_res2, x = par_opt, method.args = list(zero.tol = .Machine$double.eps, r=6))
  }
  list(par = get_par(theta, Temp = Temp),
       value = obj$fn(phi),
       cov = cov, jac = jac,
       exitflag = exitflag)
}
