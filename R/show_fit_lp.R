#' SHOW fit method log-periodogram (LP). Fit SHOW via LP method.
#'
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
#'   \item{`exitflag`}{The exit flag given by the optimizer, for debug purposes}
#' }
show_fit_lp <- function(fseq, Ypsd, fs, Temp,
                        bin_size, bin_type, phi0,
                        fit_type = c("direct", "incremental"),
                        optimizer = c("optim", "Adam"),
                        getHessian = FALSE, ...) {
  fit_type <- match.arg(fit_type)
  optimizer <- match.arg(optimizer)
  fbar <- binning(fseq, bin_size = bin_size, bin_type = bin_type)
  Zbar <- log(binning(Ypsd, bin_size = bin_size, bin_type = bin_type))
  constZ <- mean(Zbar) # normalize to avoid numerical overflow
  obj <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                    method = "LP_nlp",
                                    fbar = as.matrix(fbar),
                                    Zbar = as.matrix(Zbar - constZ),
                                    fs = fs/exp(bin_factor(bin_size))),
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
  theta <- c(exp(phi), tau = exp(obj$simulate(phi)$zeta + constZ))
  # numerical hessian
  he <- NULL
  if(getHessian) {
    obj_nll <- TMB::MakeADFun(data = list(model = "SHOW_nat",
                                    method = "LP_nll",
                                    fbar = as.matrix(fbar),
                                    Zbar = as.matrix(Zbar - constZ),
                                    fs = fs/exp(bin_factor(bin_size))),
                        parameters = list(phi = as.matrix(c(0,0,0)), zeta = 0),
                        silent = TRUE, DLL = "realPSD_TMBExports")
    phi_zeta <- c(exp(phi), obj$simulate(phi)$zeta)
    he <- numDeriv::hessian(func = obj_nll$fn, x = phi_zeta) 
    he <- he[1:3, 1:3] # truncate the row and col wrt tau
    # cov <- solve(he)
  }
  list(par = get_par(theta, Temp = Temp),
       value = obj$fn(phi), he = he,
       exitflag = exitflag)
}
