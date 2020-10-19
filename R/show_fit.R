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
#' @param optimizer Either "optim" (in R) or "Adam" (supplied by realPSD). For now, Adam's hyperparameter tunning is not fully supported. Learning rate and nsteps are preset. `Adam` is only used for `direct` fitting.
#' @param vcov If TRUE, return the numerical variance-covariance matrix for all the parameters (original scale)
#' @param get_jac If TRUE, return the numerical Jacobian of the SHOW parameters k, f0, Q, only valid if method == "NLS".
#' @param ... Additional arguments to [stats::optim()].
#'
#' @return A list with elements 
#' \describe{
#'   \item{`par`}{The fitted parameter value which is the minimizer of the objective function.}
#'   \item{`Temp`}{The temperature constant used in calculation, Kelvin.}
#'   \item{`value`}{The value of the objective function (negative loglikelihood) at the given `par`.}
#'   \item{`cov`}{The numerical variance-covariance matrix for all the parameters (original scale).}
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
                     optimizer = c("optim", "Adam"),
                     vcov = FALSE,
                     get_jac = FALSE, ...) {
  Ypsd <- Ypsd * fs # since Ypsd is returned by periodogram which has been scaled by 1/fs.
  method <- match.arg(method)
  if(method == "lp") {
    out <- show_fit_lp(fseq = fseq, Ypsd = Ypsd, fs = fs, Temp = Temp,
                       bin_size = bin_size, bin_type = bin_type,
                       phi0 = phi0, fit_type = fit_type, optimizer = optimizer,
                       vcov = vcov, ...)
  } else if(method == "nls") {
    out <- show_fit_nls(fseq = fseq, Ypsd = Ypsd, fs = fs, Temp = Temp,
                        bin_size = bin_size, bin_type = bin_type,
                        phi0 = phi0, fit_type = fit_type, optimizer = optimizer,
                        vcov = vcov, get_jac = get_jac, ...)
  } else if(method == "mle") {
    out <- show_fit_mle(fseq = fseq, Ypsd = Ypsd, fs = fs, Temp = Temp,
                        phi0 = phi0, fit_type = fit_type, optimizer = optimizer,
                        vcov = vcov, ...)
  }
  out
}
