# check sim_time using OU process.

#' Generate OU process at regular time intervals.
#'
#' @param gamma,mu,sigma OU process parameters (scalars).
#' @param dT Interobservation time (scalar).
#' @param x0 Initial values.  If missing sampled from the OU stationary distribution.
#' @param nsteps Number of steps to simulate the process.
#' @param nts Number of iid trajectories to generate.
#' @return `nsteps x nts` matrix of OU realizations.
#' @details The OU process is defined as
#' ```
#' dx_t = -gamma x_t dt + sqrt(2*gamma) dB_t.
#' ```
#' It has mean `E[x_t] = mu`, autocorrelation `cor(x_0, x_t) = e^{-gamma*t}`, and long-run variance `var(x_t) = sigma^2/(2*gamma)`.
#'
#' For equally-spaced time points, the OU process can be generated from an AR(1) process, which can be efficiently calculated without for-loops using the built-in function [stats::filter()].
ou_sim <- function(gamma, mu, sigma, dT, x0, nsteps, nts = 1) {
  # initial values
  if(missing(x0)) {
    x0 <- rnorm(nts, mean = mu, sd = sigma/sqrt(2*gamma))
  } else {
    x0 <- rep(x0, len = nts)
  }
  # standard deviation of AR noise
  osd <- sigma*sqrt((1-exp(-2*gamma*dT))/(2*gamma))
  ar_coef <- exp(-gamma*dT) # AR filter coefficient
  ## eps <- ts(matrix(rnorm(nsteps*nts, sd = sd), nsteps, nts),
  ##                  start = 0, dT = dT) # used ts
  eps <- matrix(rnorm(nsteps*nts, sd = osd), nsteps, nts)
  ## if(with_ts) eps <- ts(eps, start = 0, deltat = dT)
  ## eps[1,] <- x0
  eps <- filter(x = eps, filter = ar_coef, init = t(x0),
                       method = "recursive")
  eps <- eps + c(0, filter(rep(mu*(1-ar_coef), nsteps-1),
                           ar_coef, "recursive"))
  tsp(eps) <- NULL # remove time series attributes
  if(nts == 1) eps <- c(eps) # return vector for single series
  eps
}


# example: OU process.
# dX_t = -gamma X_t dt + sqrt(2*gamma) dB_t
# acf(t) = exp(-gamma*|t|)

# continuous time
ou_cpsd <- function(fseq, gamma) {
  1/(gamma^2 + 4*(pi*fseq)^2)
}

## # discrete time approximation to SDE
## # X_n+1 = (1 - gamma*dt) * X_n + sqrt(2*gamma*dT) * Z_n,  Z_n ~iid N(0,1)
## # => acf_X(n) = 2*gamma*dT/(1-phi^2) * phi^|n|,  phi = (1-gamma*dt)
## ou_dpsd <- function(fseq, gamma, fs) {
##   omega <- 2*pi*fseq/fs
##   phi <- (1-gamma/fs)
##   2*gamma/fs/(1+phi^2 - 2*phi*cos(omega))
## }

ou_dpsd <- function(fseq, gamma, fs) {
  phi <- exp(-gamma/fs)
  pcv <- 2*phi * cos(2*pi * fseq/fs)
  ## nu <- 2*pi * fseq/fs
  ## p1 <- 1/(1-exp(-1i*nu - gamma/fs)) + 1/(1-exp(1i*nu - gamma/fs)) - 1
  ## p2 <- (2 - pcv)/(phi^2 - pcv + 1) - 1
  (1-phi^2)/(phi^2 - pcv + 1)/(2*gamma)
  ## (2 - pcv)/(phi^2 - pcv + 1)/(2*gamma)
  ## (1-phi^2)/(phi^2 - pcv + 1)/(2*gamma)
}

# check how these match up
gamma <- 10 # true gamma
# frequency range of interest
fmin <- 5e-2
fmax <- 5e+3
# sampling frequency for discrete-time _approximation_
fsamp <- 2 * fmax # Hz

curve(fsamp * ou_cpsd(fseq = x, gamma = gamma),
      from = fmin, to = fmax, log = "xy",
      xlab = "Frequency (Hz)",
      ylab = "PSD")
curve(ou_dpsd(fseq = x, gamma = gamma, fs = fsamp),
      from = fmin, to = fmax, add = TRUE, col = "red")

# check variance
2 * integrate(ou_cpsd, lower = 0, upper = Inf, gamma = gamma)$val
1/(2*gamma)

# ok.  simulate data, plot fft

gamma <- exp(runif(1, -3, 2))
sigma <- exp(runif(1, -3, 2))
mu <- 0
fs <- 1e2 * (1 + runif(1, -.05, .05))
N <- 1e6
## sim_domain <- "time"
sim_domain <- "freq"
if(sim_domain == "time") {
  xt <- ou_sim(gamma = gamma, mu = mu, sigma = sigma, dT = 1/fs,
               nsteps = N)
} else if(sim_domain == "freq") {
  fseq <- seq(0, fs/2, len = 1e4)
  xt <- sim_time(fs = fs, N = N, fseq = fseq,
                 cpsd = sigma^2 * ou_cpsd(fseq = fseq, gamma = gamma))
}

# calculate psd
Yt <- abs(fftw::FFT(xt)[1+1:(N/2)])^2/N
fseq <- 1:(N/2)/N * fs
bin_size <- 1000
fbin <- binning(x = fseq, bin_size = bin_size)
epsd <- binning(x = Yt, bin_size = bin_size)

plot(fbin, epsd, type = "l", log = "xy")
lines(fbin, sigma^2 * fs * ou_cpsd(fseq = fbin, gamma = gamma),
      col = "red")
lines(fbin, sigma^2 * ou_dpsd(fseq = fbin, gamma = gamma, fs = fs),
      col = "blue")


#--- ok just with white noise --------------------------------------------------

sigma <- rexp(1)
N <- 5e6
xt <- rnorm(N, sd = sigma)

Yt <- abs(fftw::FFT(xt)[1:(N/2)])^2/N
bin_size <- 10000
epsd <- binning(x = Yt, bin_size = bin_size)

plot(epsd, type = "l", log = "xy")
abline(h = sigma^2, col = "red")
