#' Calculate SHOWF PSD.
#'
#' @param fseq Frequency vector (units: Hz).
#' @param k Cantilever stiffness (N/m).
#' @param f0 Resonance frequency (Hz).
#' @param Q Quality factor (unitless).
#' @param Sw White noise power (m^2/Hz).
#'
#' @return The SHOW PSD evaluated at `fseq` (m^2/Hz).
#' @note Basically the same as psdSHO.R but psdSHO offers the choice of unit conversion from m2/Hz to fm2/Hz.
#' @export
showf_psd <- function(fseq, k, f0, Q, Sw, Af, alpha, Temp) {
  # theta <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
  ## Kb <- 1.381e-23             # Boltzmann's constant
  ## tau <- Kb*Temp/(k*pi*f0*Q)
  ## Rw <- Sw/tau
  # f2 <- (fseq/theta["f0"])^2 
  # theta["tau"] * (theta["Rw"] + 1/((1 - f2)^2 + f2/theta["Q"]^2))
  show_psd(fseq, k, f0, Q, Sw, Temp) + Af/fseq^alpha
}
