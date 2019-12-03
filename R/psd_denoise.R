#' Remove periodic noise from the periodogram.
#'
#' @param Yf Vector of periodogram ordinates.
#' @param psd_fs Discrete-time PSD at the values in `Yf`, i.e., theoretical expectation of `Yf` in the absence of noise.
#' @param alpha Type-I error tolerance level.
#'
#' @return The vector `Yf` with outliers replaced by random Exponential draws with mean given by `psd_d`.
#' @export
psd_denoise <- function(Yf, psd_fs, alpha = .01) {
  nY <- length(Yf) # number of periodogram ordinates
  Zf <- Yf/psd_fs # normalized periodogram ordinates
  # do while loop (implemented as for-loop + break)
  iout <- NULL # indices of outliers
  for(ii in 1:nY) {
    imax <- which.max(Zf)
    Gstat <- Zf[imax]/sum(Zf) # observed value of G
    pval <- fisherG(a = Gstat, q = nY)
    if(pval < alpha) {
      # found outlier
      iout <- c(iout, imax)
      Zf[imax] <- rexp(1) # replace by a random exponential
    } else break
  }
  # replace periodogram values
  if(length(iout) > 0) Yf[iout] <- psd_fs[iout] * Zf[iout]
  Yf
}

#--- depreciated ---------------------------------------------------------------

## psd_denoise <- function(fseq, psd_noise, Q, f0, tau, Rw, freq_range) {
##   # determine range of data to fit
##   f_lb <- freq_range[1]
##   f_ub <- freq_range[2]
##   cond <- which(fseq >= f_lb & fseq <= f_ub)
##   fseq_cond <- fseq[cond]
##   psd_cond <- psd_noise[cond]
##   # remove periodic components from PSD
##   q <- length(psd_cond)
##   acut <- 1/q
##   psd_fit <- exp(logSHOW(fseq_cond, f0, Q, tau, Rw))
##   gpsd <- psd_cond / psd_fit
##   a <- max(gpsd/sum(gpsd))
##   ind <- which(gpsd/sum(gpsd) == a)
##   g <- fisherGstat(a, q)
##   ncut <- 0
##   while(g < acut) {
##     ncut <- ncut + 1
##     psd_cond[ind] <- rexp(length(ind),
##       rate = 1/exp(logSHOW(fseq_cond[ind], f0, Q, tau, Rw)))
##     gpsd <- psd_cond / psd_fit
##     a <- max(gpsd/sum(gpsd))
##     ind <- which(gpsd/sum(gpsd) == a)
##     g <- fisherGstat(a, q)
##   }
##   yPSD_clean <- psd_noise
##   yPSD_clean[cond] <- psd_cond
##   return(yPSD_clean)
## }
