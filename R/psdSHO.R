#' @title Power spectral density of SHO model.
#' @param f Frequency or a vector of frequencies in Hz.
#' @param f_0 Resonance frequency f_0.
#' @param Q Quality factor.
#' @param k Cantilever stiffness.
#' @param Kb Boltzmann constant.
#' @param T Temperature in Kelvin.
#' @param unit_conversion Indicator of whether or not to convert the unit of PSD to fm2/Hz. The femtometer (fm) is an SI unit of length equal to 1e-15 meters.
#' @export
psdSHO <- function(f, f0, Q, k, Kb, T, unit_conversion = TRUE) {
  numerator <- Kb * T / (k * pi * f0 * Q)
  denominator <- ((f/f0)^2 - 1)^2 + (f/(f0*Q))^2
  psd <- numerator/denominator
  CONSTANT <- 1
  if(unit_conversion == TRUE) CONSTANT = 1e30
  return(psd * CONSTANT)
}
