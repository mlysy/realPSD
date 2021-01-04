#' Calculate SHOWF PSD.
#'
#' @param freq Frequency vector (units: Hz).
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw Power of the white noise term (m^2/Hz).
#' @param Sf Power of the 1/f noise term (m^2/Hz).
#' @param alpha Exponent of the 1/f noise term (unitless).
#' @param Temp Temperature (Kelvin).
#'
#' @return The SHOWF PSD evaluated at `freq` (m^2/Hz).
#' @export
showf_psd <- function(freq, k, f0, Q, Sw, Sf, alpha, Temp) {
  # theta <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
  ## Kb <- 1.381e-23             # Boltzmann's constant
  ## tau <- Kb*Temp/(k*pi*f0*Q)
  ## Rw <- Sw/tau
  # f2 <- (freq/theta["f0"])^2
  # theta["tau"] * (theta["Rw"] + 1/((1 - f2)^2 + f2/theta["Q"]^2))
  show_psd(freq, k, f0, Q, Sw, Temp) + Sf/freq^alpha
}
