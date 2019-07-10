#' @title Unwrap radian phases by adding multiples of 2*pi as appropriate to
#' remove jumps greater than tol. tol defaults to pi.
#' @param a Radian angle.
#' @param tol Tolerance.
#' @param dim Dimension.
#' @export
unwrap <- function(a, tol = pi, dim = 1) {
  sz <- dim(a)
  nd <- length(sz)
  if (nd == 0) {
    sz <- length(a)
    nd <- 1
  }
  if (! (length(dim) == 1 && dim == round(dim)) && dim > 0 && dim < (nd + 1))
    stop("unwrap: dim must be an integer and valid dimension")
  # Find the first non-singleton dimension
  while (dim < (nd + 1) && sz[dim] == 1)
    dim <- dim + 1
  if (dim > nd)
    dim <- 1
  # Don't let anyone use a negative value for TOL.
  tol <- abs(tol)
  rng <- 2*pi
  m <- sz[dim]
  # Handle case where we are trying to unwrap a scalar, or only have
  # one sample in the specified dimension.
  if (m == 1)       
    return(a)
  # Take first order difference to see so that wraps will show up
  # as large values, and the sign will show direction.
  idx <- list()
  for (i  in  1:nd) 
    idx[[i]] <- 1:sz[i]
  idx[[dim]] <- c(1,1:(m-1))
  d <- a[unlist(idx)] - a
  # Find only the peaks, and multiply them by the range so that there
  # are kronecker deltas at each wrap point multiplied by the range
  # value.
  p <-  rng * (((d > tol) > 0) - ((d < -tol) > 0))
  # Now need to "integrate" this so that the deltas become steps.
  if (nd == 1)
    r <- cumsum(p)
  else  
    r <- apply(p, MARGIN = dim, FUN = cumsum)
  # Now add the "steps" to the original data and put output in the
  # same shape as originally.
  a + r
} 

# test P1 == P2 == P3
# Z = c(1 - 1i, 2 + 1i, 3 - 1i, 4 + 1i, 
#   1 + 2i, 2 - 2i, 3 + 2i, 4 - 2i, 
#   1 - 3i, 2 + 3i, 3 - 3i, 4 + 3i, 
#   1 + 4i, 2 - 4i, 3 + 4i, 4 - 4i)
# P1 = Arg(Z)
# P2 = atan2(Im(Z), Re(Z))
# P3 = Im(log(Z))

#' @title noise_sin: This function returns one-sided, unfolded, FFT of sine, with DC component removed. 
#' @param A Amplitude.
#' @param f Frequency.
#' @param T Total time.
#' @param SF Sampling Freq (Hz).
#' @param detrend_flag 1 = True, 0 = False.
#' @return A list of 
#' xTime: Time values, 
#' yTime: Sine signal values, 
#' xFreq: Frequencies,
#' yFreq: First freq. up to Nyquist freq of: fft(yTime)/(sqrt(FR)*N)
#' @details NOTE: Removes DC component, select first frequency up to the Nyquist frequency. NOTE: Uses even PSD's only so it will round yTime to multiple of 2
noise_sin <- function(f_noise, A_noise, T, SF, detrend_flag) {
  ## Begin
  # Generate sine-wave in time domain
  # xTime = [1/SF:1/SF:T]
  xTime <- seq(from = 1/SF, to = T, by = 1/SF)
  yTime <- A_noise * sin(2*pi*xTime*f_noise)

  # FFT to frequency domain
  if(detrend_flag == 1)
      yTime <- detrend(yTime) # detrend
  N <- length(yTime) - length(yTime)%%2 # use even PSD's only
  yTime <- yTime[1:N]
  FR <- 1/T # frequency resolution (Hz)

  yFreq <- fftw::FFT(yTime)
  # remove DC component, select first frequency up to the Nyquist frequency.
  # drop negative frequencies. For a real signal, negative frequencies are redundant.
  yFreq <- yFreq[2: (N/2+1)] # underwent folding so must double to counter loss of power
  xFreq <- seq(from = FR, to = SF/2, length.out = N/2)
  # assert(length(xFreq) == length(yFreq))
  if(length(xFreq != yFreq)) stop("The length of xFreq must be equal to that of yFreq!")

  return(list(xTime = xTime, yTime = yTime, xFreq = xFreq, yFreq = yFreq))
}
# log10.axis <- function(side, at, ...) {
#     at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
#     lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
#     axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
#     axis(side=side, at=at, labels=lab, ...)
# }
# minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){

#   lims <- par("usr")
#   if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

#   major.ticks <- pretty(lims,n=5)
#   if(missing(mn)) mn <- min(major.ticks)
#   if(missing(mx)) mx <- max(major.ticks)

#   major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]

#   labels <- sapply(major.ticks,function(i)
#             as.expression(bquote(10^ .(i)))
#           )
#   axis(ax,at=major.ticks,labels=labels,...)

#   n <- n+2
#   minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
#   minors <- minors[-c(1,n)]

#   minor.ticks = c(outer(minors,major.ticks,`+`))
#   minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


#   axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
# }
