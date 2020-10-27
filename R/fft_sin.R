#' @title FFT (discrete Fourier transform) of sine wave noise
#' @param fseq A sequence of frequency points (must be points in the frequency basis)
#' @param f0 Resonance frequency
#' @param Q Quality factor
#' @param fs Sampling frequency
#' @param unit_conversion TRUE/FALSE, if TRUE, use fm2/Hz
#' @export
fft_sin <- function(fseq, f0, Q, fs, unit_conversion) {
  # if(is.null(N))
  N <- length(fseq)
  Const <- 1e30                # unit conversion, 1 m2 = 1e30 fm2
  if(!unit_conversion) {
    D <- (Q^0.5) * 3.5e3 / Const
  } else {
    D <- (Q^0.5) * 3.5e3 # 3.5e3 is copied from the MATLAB code
  }
  # k  <- 0.172                  # Cantilever stiffness, N/m
  # Kb <- 1.381e-23              # Boltzmann's constant, (m2*kg)/(s2*K)
  # Temp <- 298                  # Temperature, K
  # D <- max(psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion))
  dT <- 1/fs
  xi <- rnorm(1, f0, 10)
  phi <- runif(1, 0, 2*pi)
  sin_fft <- D/(2*1i * sqrt(N)) * (
    exp(phi*1i) * (exp(2*pi*1i*(xi-fseq)*dT*N)-1)/(exp(2*pi*1i*(xi-fseq)*dT)-1) -
    exp(-phi*1i) * (exp(-2*pi*1i*(xi+fseq)*dT*N)-1)/(exp(-2*pi*1i*(xi+fseq)*dT)-1)
  )
  return(sin_fft)
}
