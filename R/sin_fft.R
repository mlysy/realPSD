#' @title FFT (discrete Fourier transform) of sine wave noise
#' @param fseq A sequence of frequency points (must be points in the frequency basis)
#' @param fs Sampling frequency
#' @param D Amplitude
#' @param xi as specified in the paper 
#' @param phi as specified in the paper
#' @export
sin_fft <- function(fseq, fs, D, xi, phi) {
  # if(is.null(N))
  N <- length(fseq)
  dT <- 1/fs
  sin_fft <- D/(2*1i) * (
    exp(phi*1i) * (exp(2*pi*1i*(xi-fseq)*dT*N)-1)/(exp(2*pi*1i*(xi-fseq)*dT)-1) -
    exp(-phi*1i) * (exp(-2*pi*1i*(xi+fseq)*dT*N)-1)/(exp(-2*pi*1i*(xi+fseq)*dT)-1)
  )
  return(sin_fft)
}
