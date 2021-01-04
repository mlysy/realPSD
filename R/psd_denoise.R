#' Remove periodic noise from the periodogram.
#'
#' @param Ypsd Vector of periodogram ordinates.
#' @param psd PSD at the values in `Ypsd`, i.e., theoretical expectation of `Ypsd` in the absence of noise.
#' @param alpha Type-I error tolerance level.
#' @param method Denoising method to use.  Either `fisherG` for Fisher's G-statistic test, or `BH` for the Benjamini-Hochberg method, in which case `alpha` specifies the false discovery rate.
#'
#' @return The vector `Ypsd` with outliers replaced by random Exponential draws with mean given by `psd`.
#' @export
psd_denoise <- function(Ypsd, psd, alpha = .01,
                        method = c("fisherG", "BH")) {
  method <- match.arg(method)
  ## nY <- length(Ypsd) # number of periodogram ordinates
  Zf <- Ypsd/psd # normalized periodogram ordinates
  ## # do while loop (implemented as for-loop + break)
  ## iout <- NULL # indices of outliers
  ## for(ii in 1:nY) {
  ##   imax <- which.max(Zf)
  ##   Gstat <- Zf[imax]/sum(Zf) # observed value of G
  ##   pval <- fisherG(a = Gstat, q = nY)
  ##   if(pval < alpha) {
  ##     # found outlier
  ##     iout <- c(iout, imax)
  ##     Zf[imax] <- rexp(1) # replace by a random exponential
  ##   } else break
  ## }
  if(method == "fisherG") {
    dn <- fisherG_denoise(Zf, alpha = alpha)
  } else if(method == "BH") {
    dn <- BH_denoise(Zf, alpha)
  }
  # replace periodogram values
  ## if(length(iout) > 0) Ypsd[iout] <- psd[iout] * Zf[iout]
  if(!is.null(dn)) {
    Ypsd[dn$iout] <- psd[dn$iout] * dn$Zf
  }
  Ypsd
}

#--- helper functions ----------------------------------------------------------

#' Denoising with Fisher-G method.
#'
#' @noRd
fisherG_denoise <- function(Zf, alpha) {
  nZ <- length(Zf) # number of periodogram ordinates
  Grange <- fisherG_range(alpha, Gstat = max(Zf)/sum(Zf) / 2^(0:20), q = nZ)
  root <- uniroot(f = function(logG) log(fisherG(exp(logG), nZ)) - log(alpha),
                  interval = log(Grange))$root
  Gcut <- exp(root)
  iout <- Zf/sum(Zf) > Gcut
  list(iout = iout, Zf = rexp(sum(iout)))
  ## ## fisherG(a = Gstat, q = nZ)
  ## ## Gstat * exp(-2:2)
  ## Gpval <- sapply(Gstat / scale, fisherG, q = nZ)
  ## nG <- length(Gpval)
  ## ikeep <- !is.na(Gpval)
  ## ikeep[ikeep] <- Gpval[ikeep] > 0
  ## ikeep[ikeep] <- Gpval[ikeep] <= 1
  ## keep1 <- 1:nG < which.max(ikeep) # up to first keep
  ## # remove after first fail
  ## ikeep[which.min(keep1 | ikeep) < 1:nG] <-  FALSE
  ## ikeep[keep1] <- FALSE # remove before first keep
  ## if(!any(ikeep)) stop("Invalid scale.")
  ## ## Gstat <- max(Zf/sum(Zf))
  ## ## ikeep <- maxbad-(2:1)
  ## ## ikeep <- 1:nZ
  ## ## ikeep[1:(maxbad-1)] <- NA
  ## Gdata <- data.frame(scale = log2(scale[ikeep]), pval = log(Gpval[ikeep]))
  ## # take slope of last two points
  ## i <- nrow(Gdata)-1
  ## slope <- (Gdata$pval[i+1]-Gdata$pval[i])/(Gdata$scale[i+1]-Gdata$scale[i])
  ## fisherG(Gstat/exp(Gdata$scale[i+1] - Gdata$pval[i+1]/slope), q = nZ)
  ## plot(Gdata$scale, Gdata$pval)
  ## lines(Gdata$scale,
  ##       predict(lm(pval ~ scale + I(scale^2), data = Gdata)), col = "red")
  ## sapply((Gstat * scale)[1:(maxbad-1)], fisherG, q = nZ)
  ## plot()
}

## fisherG_denoise <- function(Zf, alpha) {
##   nZ <- length(Zf) # number of periodogram ordinates
##   # do while loop (implemented as for-loop + break)
##   iout <- NULL # indices of outliers
##   for(ii in 1:nZ) {
##     imax <- which.max(Zf)
##     Gstat <- Zf[imax]/sum(Zf) # observed value of G
##     pval <- fisherG(a = Gstat, q = nZ)
##     if(pval < alpha) {
##       # found outlier
##       iout <- c(iout, imax)
##       Zf[imax] <- rexp(1) # replace by a random exponential
##     } else break
##   }
##   if(is.null(iout)) return(NULL)
##   list(iout = iout, Zf = Zf[iout])
## }

#' Denoising with Benjamini-Hochberg method.
#'
#' @noRd
BH_denoise <- function(Zf, alpha) {
  nZ <- length(Zf) # number of periodogram ordinates
  # do procedure exactly
  Uk <-  1 - exp(-Zf)
  iord <- order(Uk)
  Uk <- Uk[iord]
  kstar <- 1 - Uk <= alpha * (nZ - 1:nZ + 1)/nZ
  if(!any(kstar)) return(NULL) # no outliers
  # all Uk greater than or equal to this rank get replaced
  kstar <- which.max(kstar)
  Umin <- ifelse(kstar > 1, Uk[kstar-1], 0)
  nout <- nZ - kstar + 1
  Ustar <- sort(runif(n = nout, min = Umin, max = 1))
  Zstar <- -log(1-Ustar)
  list(iout = tail(iord, nout), Zf = Zstar)
}

fisherG_range <- function(alpha, Gstat, q) {
  Gstat <- sort(Gstat)
  Gpval <- sapply(Gstat, fisherG, q)
  ikeep <- !is.na(Gpval) # remove NAs
  # keep only values between 0 and 1
  ikeep[ikeep] <- Gpval[ikeep] > 0
  ikeep[ikeep] <- Gpval[ikeep] <= 1
  if(!any(Gpval[ikeep] < alpha) || !any(Gpval[ikeep] > alpha)) {
    stop("Gstat does not give p-values below and above alpha.")
  }
  c(Gstat[ikeep][which.max(alpha < Gpval[ikeep])],
    Gstat[ikeep][which.max(Gpval[ikeep] < alpha)])
}

#--- scratch -------------------------------------------------------------------

## fbar <- binning(psd_fit$freq, 100)
## Ybar <- binning(psd_fit$Ypsd, 100)
## Ybar2 <- binning(Ypsd, 100)
## pbar <- binning(psd, 100)
## plot(fbar, Ybar, type = "l", log = "xy")
## lines(fbar, pbar, col = "red")

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
