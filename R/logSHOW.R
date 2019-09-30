#' @title log Power spectral density of continuous-time SHO model.
#' @param fseq Frequency or a vector of frequencies in Hz.
#' @param f_0 Resonance frequency f_0.
#' @param Q Quality factor.
#' @param tau kB * Temp / (k*pi*f0*Q)
#' @param Rw Aw/tau
#' @param unit_conversion Indicator of whether or not to convert the unit of PSD to fm2/Hz. The femtometer (fm) is an SI unit of length equal to 1e-15 meters.
#' @export
logSHO <- function(fseq, f0, Q, tau, Rw, unit_conversion = TRUE) {
  log_psd <- log(tau) + log(Rw + 1/((1 - (fseq/f0)^2)^2 + (fseq/(f0*Q))^2))
  CONSTANT <- 1
  if(unit_conversion == TRUE) CONSTANT = 1e30
  log_psd <- log_psd + log(CONSTANT)
  return(log_psd)
}
