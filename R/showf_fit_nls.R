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
#' @param getHessian If TRUE, return the numerical Hessian matrix of the original SHOW parameter f0, Q, Rw
#' @param ... Additional arguments to [stats::optim()].
#' 
#' @return A list with elements 
#' \describe{
#'   \item{`par`}{The fitted parameter value which is the minimizer of the objective function.}
#'   \item{`Temp`}{The temperature constant used in calculation, Kelvin.}
#'   \item{`value`}{The value of the objective function (negative loglikelihood) at the given `par`.}
#'   \item{`he`}{The numerical Hessian matrix of the SHOW model parameter f0, Q, Rw.}
#'   \item{`jac`}{The numerical Jacobian matrix of the SHOW model parameter f0, Q, Rw.}
#'   \item{`exitflag`}{The exit flag given by the optimizer, for debug purposes}
#' }
showf_fit_nls <- function(fseq, Ypsd, fs, Temp,
                         bin_size, bin_type, phi0,
                        fit_type = c("direct", "incremental"),
                        optimizer = c("optim", "Adam"),
                        getHessian = FALSE, getJacobian = FALSE, ...) {
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
  # hessian
  he <- NULL
  if(getHessian) {    
    obj_nll <- TMB::MakeADFun(data = list(model = "SHOWF_nat",
                                    method = "NLS_nll",
                                    fbar = as.matrix(fseq),
                                    Ybar = as.matrix(Ypsd/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(rep(0,5)), tau = 0),
                        map = map,
                        silent = TRUE, DLL = "realPSD_TMBExports")
    phi_tau <- c(exp(phi), tau = obj$simulate(phi)$tau)
    he <- numDeriv::hessian(func = obj_nll$fn, x = phi_tau)
    # browser()
    # ind <- c(1,2,4,5)
    # tmp <- he[ind, ind]
    # solve(tmp)
    he <- he[c(1,2,5), c(1,2,5)] # truncate the row and col wrt tau
    # cov <- solve(he)
  }
  # Jacobian
  jac <- NULL
  if(getJacobian) {
    obj_res <- TMB::MakeADFun(data = list(model = "SHOWF_nat",
                                    method = "NLS_res",
                                    fbar = as.matrix(fbar),
                                    Ybar = as.matrix(Ybar/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(rep(0,5))),
                        map = map,
                        silent = TRUE,
                        ADreport = TRUE,
                        DLL = "realPSD_TMBExports")
    nls_res2 <- function(phi) setNames((obj_res$fn(phi))^2, nm = NULL)
    phi_tau <- c(exp(phi), tau = obj$simulate(phi)$tau)
    jac <- numDeriv::jacobian(nls_res2, x = phi_tau[1:5]) # return the full Jacobian for all params in phi
    # jac <- numDeriv::jacobian(nls_res2, x = phi_tau)
  }
  list(par = get_par_showf(theta, Temp = Temp),
       value = obj$fn(phi),
       he = he, jac = jac,
       exitflag = exitflag)
}
