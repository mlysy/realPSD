#' @title PSD denoise
#' @param fseq Frequency vector
#' @param psd_noise Noisy PSD
#' @param Q Quality factor
#' @param f0 Resonance frequency
#' @param k k
#' @param Temp Temperature
#' @param unit_conversion indicator if unit conversion is required
#' @param Aw White noise amplitutde
#' @param freq_range frequency range, i.e. f0 +- f0/sqrt(2)
#' @export
psd_denoise <- function(fseq, psd_noise, Q, f0, k, Temp, unit_conversion, Aw, freq_range) {
  # determine range of data to fit
  f_lb <- freq_range[1]
  f_ub <- freq_range[2]
  cond <- which(fseq >= f_lb & fseq <= f_ub)
  fseq_cond <- fseq[cond]
  psd_cond <- psd_noise[cond]
  # remove periodic components from PSD
  q <- length(psd_cond)
  acut <- 1/q
  psd_fit <- psdSHO(fseq_cond, f0, Q, k, Temp, unit_conversion) + Aw
  gpsd <- psd_cond / psd_fit
  a <- max(gpsd/sum(gpsd))
  ind <- which(gpsd/sum(gpsd) == a)
  g <- fisherGstat(a, q)
  ncut <- 0
  while(g < acut) {
    ncut <- ncut + 1
    psd_cond[ind] <- rexp(length(ind), 
      rate = 1/(psdSHO(fseq_cond[ind], f0, Q, k, Temp, unit_conversion) + Aw))
    gpsd <- psd_cond / psd_fit
    a <- max(gpsd/sum(gpsd))
    ind <- which(gpsd/sum(gpsd) == a)
    g <- fisherGstat(a, q)
  }
  yPSD_clean <- psd_noise
  yPSD_clean[cond] <- psd_cond
  return(yPSD_clean)
}
