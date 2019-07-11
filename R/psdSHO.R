#' @title Power spectral density of SHO model.
#' @param f Frequency or a vector of frequencies in Hz.
#' @param f_0 Resonance frequency f_0.
#' @param Q Quality factor.
#' @param k Cantilever stiffness.
#' @param kb Boltzmann constant.
#' @param T Temperature in Kelvin.
#' @export
psd_sho <- function(f, f_0, Q, k, kb, T) {
  numerator <- kb * T / (k * pi * f_0 * Q)
  denominator <- ((f/f_0)^2 - 1)^2 + (f/(f_0*Q))^2
  return(numerator/denominator)
}
