#' SHOW fit method maximum likelihood estimation (MLE)
#'
#' @param fseq Frequency vector (units: Hz).
#' @param Ypsd PSD vector.
#' @param fs Sampling frequency.
#' @param Temp True temperature.
#' @param phi0 Parameter value (transformed scale) to initialize optimization.
#' @param bin_type Either "mean" or "median".  The former is more efficient, the latter is more robust, so better for preliminary SHOW fit to denoise against.
#' @param fit_type If "direct", fit all parameters at once.  If "incremental", do one, then two, then three, etc.
#' @param getHessian If TRUE, return the numerical Hessian matrix of the original SHOW parameter f0, Q, Rw
#' @param ... Additional arguments to [stats::optim()].
#' 
#' @return A list with elements 
#' \describe{
#'   \item{`par`}{The fitted parameter value which is the minimizer of the objective function.}
#'   \item{`Temp`}{The temperature constant used in calculation, Kelvin.}
#'   \item{`value`}{The value of the objective function (negative loglikelihood) at the given `par`.}
#'   \item{`he`}{The numerical Hessian matrix of the SHOW model parameter f0, Q, Rw.}
#'   \item{`exitflag`}{The exit flag given by the optimizer, for debug purposes}
#' }
show_fit_mle <- function(fseq, Ypsd, fs, Temp, phi0,
                        fit_type = c("direct", "incremental"),
                        getHessian = FALSE, ...) {
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
    # fit <- optim(par = phi,
    #              fn = obj$fn, gr = obj$gr, method = "BFGS", ...)
    fit <- adam(theta0 = phi, fn = obj$fn, gr = obj$gr, nsteps = 500, alpha = 1e-4)
    phi <- fit$par
    exitflag <- c(exitflag, fit$convergence)
  }
  # construct final estimator
  theta <- c(exp(phi), tau = obj$simulate(phi)$tau * constY)
  # cov
  # cov <- NULL
  he <- NULL
  if(getHessian) {    
    obj_nll <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                    method = "MLE_nll",
                                    f = as.matrix(fseq),
                                    Y = as.matrix(Ypsd/constY),
                                    fs = fs),
                        parameters = list(phi = as.matrix(c(0,0,0)), tau = 0),
                        silent = TRUE, DLL = "realPSD_TMBExports")
    phi_tau <- get_phi(par = get_par(theta, Temp), 
      method = "MLE", Temp = Temp, const = constY)
    # phi_tau[4] <- phi_tau[4] / fs
    he <- numDeriv::hessian(func = obj_nll$fn, x = phi_tau)
    he <- he[1:3, 1:3] # truncate the row and col wrt tau
    # cov <- solve(he)
  }
  list(par = get_par(theta, Temp = Temp),
       value = obj$fn(phi), he = he,
       exitflag = exitflag)
}