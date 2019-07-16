# to reproduce Figure 2 in the paper
require(realPSD)
# ---------- SHO model parameters ----------
T  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q_vec  <- c(1, 10, 100, 500)  # Quality factors
k  <- 0.172               # Cantilever stiffness, N/m
Kb <- 1.381e-23           # Boltzmann's constant
T <- 298                  # Temperature, Kelvin
Aw <- 19000               # white noise, fm2/Hz
# sig2 <- Kb*T/(k*pi*f0*Q)  # variance sig2
# Rw <- Aw/sig2
# alpha <- 0.55           # 1/f decay exponent
# ---------- simulate random datasets ----------
f_lb <- f0 - f0/sqrt(2) # frequency lower bound
f_ub <- f0 + f0/sqrt(2) # frequency upper bound
f <- seq(from = f_lb, to = f_ub, by = 1/T) # frequency domain, Hz
# nsim <- 1000 # number of simulations / number of simulated datasets
# system.time(
#   # tmp <- rexp(length(f)*nsim, rate = 1)
#   tmp <- rexp(length(f), rate = 1)
# )
# print(object.size(tmp), units = "GB")

# for loop
nsim <- 10
# pre-allocate space for storage
fit_Q1_lp <- matrix(NA, nsim, 3) 
fit_Q1_nls <- matrix(NA, nsim, 3) 
fit_Q1_mle <- matrix(NA, nsim, 3) 
fit_Q10_lp <- matrix(NA, nsim, 3) 
fit_Q10_nls <- matrix(NA, nsim, 3) 
fit_Q10_mle <- matrix(NA, nsim, 3) 
fit_Q100_lp <- matrix(NA, nsim, 3) 
fit_Q100_nls <- matrix(NA, nsim, 3) 
fit_Q100_mle <- matrix(NA, nsim, 3) 
fit_Q500_lp <- matrix(NA, nsim, 3) 
fit_Q500_nls <- matrix(NA, nsim, 3) 
fit_Q500_mle <- matrix(NA, nsim, 3) 
system.time(
for(ii in 1:nsim) {
  # generate a vector of exponential random variables 
  # which is the same for all Q1, Q10, Q100, Q500 but different for each iteration/simulation 
  rfreq <- rexp(length(f), rate = 1) 
  # ----- fitting simulated data -----
  # specify the Q factor
  Q1 <- Q_vec[1] 
  fit_Q1_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q1_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q1_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood
  # specify the Q factor
  Q10 <- Q_vec[2]
  fit_Q10_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q10_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q10_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood  
  # specify the Q factor
  Q100 <- Q_vec[3]
  fit_Q100_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q100_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q100_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood  
  # specify the Q factor
  Q500 <- Q_vec[4]
  fit_Q500_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q500_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q500_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood  
}
)


# nbins <- ceiling(length(f) / binSize)
# nleft <- nbins*binSize - length(f) # number of positions left
# ftmp <- c(f, rep(0, nleft)) # fill the remaining positions with 0s
# # binning the frequency domain and periodogram values
# fmat <- matrix(ftmp, nrow = binSize) # each column is a bin whose values have consecutive memory locations
# fbar <- colMeans(fmat)
# fbar[nbins] <- sum(f[(binSize*(nbins-1)+1):length(f)])/(length(f) - binSize*(nbins-1))


