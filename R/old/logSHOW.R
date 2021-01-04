#' @title log Power spectral density of continuous-time SHO model.
#' @param fseq Frequency or a vector of frequencies in Hz.
#' @param f0 Resonance frequency f0.
#' @param Q Quality factor.
#' @param tau kB * Temp / (k*pi*f0*Q), tau should supplied with unit_conversion if it is required
#' @param Rw Aw/tau
#' @export
logSHOW <- function(fseq, f0, Q, tau, Rw) {
  log_psd <- log(tau) + log(Rw + 1/((1 - (fseq/f0)^2)^2 + (fseq/(f0*Q))^2))
  return(log_psd)
}
