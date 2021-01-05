#' Calculate SHOW PSD.
#'
#' @param freq Frequency vector (units: Hz).
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw White noise power (m^2/Hz).
#' @param Temp Temperature (Kelvin).
#'
#' @return The SHOW PSD evaluated at `freq` (m^2/Hz).
#' @note Basically the same as psdSHO.R but psdSHO offers the choice of unit conversion from m2/Hz to fm2/Hz.
#' @export
show_psd <- function(freq, k, f0, Q, Sw, Temp) {
  theta <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
  ## Kb <- 1.381e-23             # Boltzmann's constant
  ## tau <- Kb*Temp/(k*pi*f0*Q)
  ## Rw <- Sw/tau
  f2 <- (freq/theta["f0"])^2
  theta["tau"] * (theta["Rw"] + 1/((1 - f2)^2 + f2/theta["Q"]^2))
}

#--- helper functions ----------------------------------------------------------

#' Get `theta` from SHOW parameters.
#'
#' The parameter vector `theta` is defined as
#' ```
#' theta = (f0, Q, Rw = Sw/tau, tau),
#' ```
#' where `tau = KbT/(pi k f0 Q)`.
#'
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw White noise power (m^2/Hz).
#' @param Temp Temperature (Kelvin).
#'
#' @return Numeric vector with named elements `(f0, Q, Rw, tau)`.
#' @noRd
get_theta <- function(k, f0, Q, Sw, Temp) {
  Kb <- 1.381e-23             # Boltzmann's constant
  tau <- Kb*Temp/(k*pi*f0*Q)
  Rw <- Sw/tau
  setNames(c(f0, Q, Rw, tau), nm = c("f0", "Q", "Rw", "tau"))
}
