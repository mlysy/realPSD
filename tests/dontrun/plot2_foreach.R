# doParallel and foreach version
# ---------- load required packages ----------
require(realPSD)
require(TMB)
require(tidyverse)
require(foreach) # parallel for-loop computing procedure
require(doParallel) # parallel backend for the foreach package
require(doRNG) # for reproducible foreach %dopar% loop
# ---------- SHO model parameters ----------
T  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q_vec  <- c(1, 10, 100, 500)  # Quality factors
k  <- 0.172               # Cantilever stiffness, N/m
T <- 298                  # Temperature, Kelvin
Aw <- 19000               # white noise, fm2/Hz
# ---------- generate frequency domain ----------
f_lb <- f0 - f0/sqrt(2) # frequency lower bound
f_ub <- f0 + f0/sqrt(2) # frequency upper bound
f <- seq(from = f_lb, to = f_ub, by = 1/T) # frequency domain, Hz
# ---------- simulation setup ----------
# number of simulations
nsim <- 2
# number of bin size
binSize <- 100
# specify the Q factor values
Q1 <- Q_vec[1] 
Q10 <- Q_vec[2]
Q100 <- Q_vec[3]
Q500 <- Q_vec[4]
# ---------- parallel procedure --------
# register a multi-core platform 
# ncores <- detectCores()
ncores <- 2
cl <- makeCluster(ncores) # by default the type is PSOCK which works on Windows, FORK type only works on Unix systems  
registerDoParallel(cl)  
# set seed in a non-invasive way (without changing the original %dopar% structure)
# registerDoRNG(123)
system.time(
fit <- foreach(ii = 1:nsim, 
  .combine = "rbind", 
  .packages = c("realPSD", "TMB")) %dopar% {
  # generate a vector of exponential random variables for this iteration
  rfreq <- rexp(length(f), rate = 1) 
  # ----- fitting simulated data -----
  # Q1 <- Q_vec[1] 
  fit_Q1_lp <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q1_nls <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q1_mle <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood
  names(fit_Q1_lp) <- c("f0", "Q", "k")
  names(fit_Q1_nls) <- c("f0", "Q", "k")
  names(fit_Q1_mle) <- c("f0", "Q", "k") 
  # Q10 <- Q_vec[2]
  fit_Q10_lp <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q10_nls <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q10_mle <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood  
  names(fit_Q10_lp) <- c("f0", "Q", "k")
  names(fit_Q10_nls) <- c("f0", "Q", "k")
  names(fit_Q10_mle) <- c("f0", "Q", "k")
  # Q100 <- Q_vec[3]
  fit_Q100_lp <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q100_nls <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q100_mle <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood  
  names(fit_Q100_lp) <- c("f0", "Q", "k")
  names(fit_Q100_nls) <- c("f0", "Q", "k")
  names(fit_Q100_mle) <- c("f0", "Q", "k")
  # Q500 <- Q_vec[4]
  fit_Q500_lp <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
    binSize = 100, method = "LP_nlp")  # log-periodogram
  fit_Q500_nls <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
    binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
  fit_Q500_mle <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
    binSize = 100, method = "MLE_nlp") # Whittle likelihood  
  names(fit_Q500_lp) <- c("f0", "Q", "k")
  names(fit_Q500_nls) <- c("f0", "Q", "k")
  names(fit_Q500_mle) <- c("f0", "Q", "k")  
  # ----- convert the estimated parameters to ratios -----
  # Q = 1
  fit_Q1_lp <- fit_Q1_lp %*% diag(c(1/f0, 1/Q1, 1/k))
  fit_Q1_nls <- fit_Q1_nls %*% diag(c(1/f0, 1/Q1, 1/k))
  fit_Q1_mle <- fit_Q1_mle %*% diag(c(1/f0, 1/Q1, 1/k))
  # Q = 10
  fit_Q10_lp <- fit_Q10_lp %*% diag(c(1/f0, 1/Q10, 1/k))
  fit_Q10_nls <- fit_Q10_nls %*% diag(c(1/f0, 1/Q10, 1/k))
  fit_Q10_mle <- fit_Q10_mle %*% diag(c(1/f0, 1/Q10, 1/k))
  # Q = 100
  fit_Q100_lp <- fit_Q100_lp %*% diag(c(1/f0, 1/Q100, 1/k))
  fit_Q100_nls <- fit_Q100_nls %*% diag(c(1/f0, 1/Q100, 1/k))
  fit_Q100_mle <- fit_Q100_mle %*% diag(c(1/f0, 1/Q100, 1/k))
  # Q = 500
  fit_Q500_lp <- fit_Q500_lp %*% diag(c(1/f0, 1/Q500, 1/k))
  fit_Q500_nls <- fit_Q500_nls %*% diag(c(1/f0, 1/Q500, 1/k))
  fit_Q500_mle <- fit_Q500_mle %*% diag(c(1/f0, 1/Q500, 1/k))
  # ----- combine the datasets --------
  # Q = 1
  fit_Q1_lp <- c(fit_Q1_lp, level = "Q = 1", method = "LP")
  fit_Q1_nls <- c(fit_Q1_nls, level = "Q = 1", method = "NLS")
  fit_Q1_mle <- c(fit_Q1_mle, level = "Q = 1", method = "MLE")
  # Q = 10
  fit_Q10_lp <- c(fit_Q10_lp, level = "Q = 10", method = "LP")
  fit_Q10_nls <- c(fit_Q10_nls, level = "Q = 10", method = "NLS")
  fit_Q10_mle <- c(fit_Q10_mle, level = "Q = 10", method = "MLE")
  # Q = 100
  fit_Q100_lp <- c(fit_Q100_lp, level = "Q = 100", method = "LP")
  fit_Q100_nls <- c(fit_Q100_nls, level = "Q = 100", method = "NLS")
  fit_Q100_mle <- c(fit_Q100_mle, level = "Q = 100", method = "MLE")
  # Q = 500
  fit_Q500_lp <- c(fit_Q500_lp, level = "Q = 500", method = "LP")
  fit_Q500_nls <- c(fit_Q500_nls, level = "Q = 500", method = "NLS")
  fit_Q500_mle <- c(fit_Q500_mle, level = "Q = 500", method = "MLE")
  # combine them together
  tmp <- rbind(
    fit_Q1_lp, fit_Q1_nls, fit_Q1_mle,
    fit_Q10_lp, fit_Q10_nls, fit_Q10_mle,
    fit_Q100_lp, fit_Q100_nls, fit_Q100_mle,
    fit_Q500_lp, fit_Q500_nls, fit_Q500_mle
  ) 
  return(tmp)
}
)
# close multi-core clusters
stopCluster(cl)
# convert the result to tibble data frame
fit <- as_tibble(fit)
# boxplot
# Q_hat / Q
ggplot(fit, aes(x = level, y = Q, fill = method)) + geom_boxplot()
# k_hat / k
ggplot(fit, aes(x = level, y = k, fill = method)) + geom_boxplot()
# f0_hat / f0
ggplot(fit, aes(x = level, y = f0, fill = method)) + geom_boxplot()

set.seed(123, kind = "L'Ecuyer-CMRG")
a <- foreach(i=1:2,.combine=cbind) %dopar% {rnorm(5)}
b <- foreach(i=1:2,.combine=cbind) %dopar% {rnorm(5)}
identical(a,b)


