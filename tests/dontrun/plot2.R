# to reproduce Figure 2 in the paper
require(realPSD)
require(TMB)
# require(tidyverse)
require(parallel)
# ---------- SHO model parameters ----------
T  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q_vec  <- c(1, 10, 100, 500)  # Quality factors
k  <- 0.172               # Cantilever stiffness, N/m
# Kb <- 1.381e-23           # Boltzmann's constant
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

# ---------- simulation ----------
nsim <- 100
binSize <- 100
# detect the number of cores
ncores <- detectCores()
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
# specify the Q factor values
Q1 <- Q_vec[1] 
Q10 <- Q_vec[2]
Q100 <- Q_vec[3]
Q500 <- Q_vec[4]
# -------- parallel version ----------
# set the seed by using L'Ecuyer-CMRG
# please see the official docs for the reason
set.seed(123, kind = "L'Ecuyer-CMRG")
# ---------- Q = 1 ---------- 
# method: log periodogram
# system.time(
fit_Q1_lp <- do.call(rbind, mclapply(1:nsim, function(ii) {
  # generate exponential random variables
  rfreq <- rexp(length(f), rate = 1)
  # save for later use
  saveRDS(rfreq, file = paste0("./data/rfreq_", ii, ".rds"))
  # fit the parameter
  # return each vector of estimated parameters as a row
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, binSize, method = "LP_nlp"), 
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
# )
saveRDS(fit_Q1_lp, file = paste0("fit_Q1_lp_", nsim, ".rds"))
# method: nonlinear least squares
# system.time(
fit_Q1_nls <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, binSize, method = "NLS_nlp"),
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
# )
saveRDS(fit_Q1_nls, file = paste0("fit_Q1_nls_", nsim, ".rds"))
# method: Whittle MLE
# system.time(
fit_Q1_mle <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, binSize, method = "MLE_nlp"),
  error = function(err)
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
# )
saveRDS(fit_Q1_mle, file = paste0("fit_Q1_mle_", nsim, ".rds"))
# ---------- Q = 10 ----------
# method: LP
fit_Q10_lp <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, binSize, method = "LP_nlp"), 
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q10_lp, file = paste0("fit_Q10_lp_", nsim, ".rds"))
# method: NLS
fit_Q10_nls <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, binSize, method = "NLS_nlp"),
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q10_nls, file = paste0("fit_Q10_nls_", nsim, ".rds"))
# method: Whittle MLE
fit_Q10_mle <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, binSize, method = "MLE_nlp"),
  error = function(err)
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q10_mle, file = paste0("fit_Q10_mle_", nsim, ".rds"))
# ---------- Q = 100 ----------
# method: LP
fit_Q100_lp <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, binSize, method = "LP_nlp"), 
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q100_lp, file = paste0("fit_Q100_lp_", nsim, ".rds"))
# method: NLS
fit_Q100_nls <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, binSize, method = "NLS_nlp"),
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q100_nls, file = paste0("fit_Q100_nls_", nsim, ".rds"))
# method: Whittle MLE
fit_Q100_mle <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, binSize, method = "MLE_nlp"),
  error = function(err)
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q100_mle, file = paste0("fit_Q100_mle_", nsim, ".rds"))
# ---------- Q = 500 -----------
# method: LP
fit_Q500_lp <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, binSize, method = "LP_nlp"), 
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q500_lp, file = paste0("fit_Q500_lp_", nsim, ".rds"))
# method: NLS
fit_Q500_nls <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, binSize, method = "NLS_nlp"),
  error = function(err) 
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q500_nls, file = paste0("fit_Q500_nls_", nsim, ".rds"))
# method: Whittle MLE
fit_Q500_mle <- do.call(rbind, mclapply(1:nsim, function(ii) {
  rfreq <- readRDS(paste0("./data/rfreq_", ii, ".rds"))
  ans <- tryCatch(
  fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, binSize, method = "MLE_nlp"),
  error = function(err)
    message(paste0("Fitting error found in the ", ii, "-th iteration."))
  )
}, mc.cores = ncores))
saveRDS(fit_Q500_mle, file = paste0("fit_Q500_mle_", nsim, ".rds"))

# for-loop sequential version
# system.time(
# for(ii in 1:nsim) {
#   # generate a vector of exponential random variables 
#   # which is the same for all Q1, Q10, Q100, Q500 but different for each iteration/simulation 
#   rfreq <- rexp(length(f), rate = 1) 
#   # saveRDS(rfreq, file = paste0("rfreq_", ii, ".rds")) # saving is slow
#   # ----- fitting simulated data -----
#   # # specify the Q factor
#   # Q1 <- Q_vec[1] 
#   fit_Q1_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
#     binSize = 100, method = "LP_nlp")  # log-periodogram
#   fit_Q1_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
#     binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
#   fit_Q1_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q1, k, T, Aw, 
#     binSize = 100, method = "MLE_nlp") # Whittle likelihood
#   # # specify the Q factor
#   # Q10 <- Q_vec[2]
#   fit_Q10_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
#     binSize = 100, method = "LP_nlp")  # log-periodogram
#   fit_Q10_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
#     binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
#   fit_Q10_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q10, k, T, Aw, 
#     binSize = 100, method = "MLE_nlp") # Whittle likelihood  
#   # # specify the Q factor
#   # Q100 <- Q_vec[3]
#   fit_Q100_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
#     binSize = 100, method = "LP_nlp")  # log-periodogram
#   fit_Q100_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
#     binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
#   fit_Q100_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q100, k, T, Aw, 
#     binSize = 100, method = "MLE_nlp") # Whittle likelihood  
#   # # specify the Q factor
#   # Q500 <- Q_vec[4]
#   fit_Q500_lp[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
#     binSize = 100, method = "LP_nlp")  # log-periodogram
#   fit_Q500_nls[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
#     binSize = 100, method = "NLS_nlp") # Nonlinear least-squares
#   fit_Q500_mle[ii, ] <- fitSHOW(f, rfreq, fs, f0, Q500, k, T, Aw, 
#     binSize = 100, method = "MLE_nlp") # Whittle likelihood  
# }
# )

# # covnert the estimated parameters to ratios
# # Q = 1
# fit_Q1_lp <- fit_Q1_lp %*% diag(c(1/f0, 1/Q1, 1/k))
# fit_Q1_nls <- fit_Q1_nls %*% diag(c(1/f0, 1/Q1, 1/k))
# fit_Q1_mle <- fit_Q1_mle %*% diag(c(1/f0, 1/Q1, 1/k))
# colnames(fit_Q1_lp) <- c("f0", "Q", "k")
# colnames(fit_Q1_nls) <- c("f0", "Q", "k")
# colnames(fit_Q1_mle) <- c("f0", "Q", "k")
# # Q = 10
# fit_Q10_lp <- fit_Q10_lp %*% diag(c(1/f0, 1/Q10, 1/k))
# fit_Q10_nls <- fit_Q10_nls %*% diag(c(1/f0, 1/Q10, 1/k))
# fit_Q10_mle <- fit_Q10_mle %*% diag(c(1/f0, 1/Q10, 1/k))
# colnames(fit_Q10_lp) <- c("f0", "Q", "k")
# colnames(fit_Q10_nls) <- c("f0", "Q", "k")
# colnames(fit_Q10_mle) <- c("f0", "Q", "k")
# # Q = 100
# fit_Q100_lp <- fit_Q100_lp %*% diag(c(1/f0, 1/Q100, 1/k))
# fit_Q100_nls <- fit_Q100_nls %*% diag(c(1/f0, 1/Q100, 1/k))
# fit_Q100_mle <- fit_Q100_mle %*% diag(c(1/f0, 1/Q100, 1/k))
# colnames(fit_Q100_lp) <- c("f0", "Q", "k")
# colnames(fit_Q100_nls) <- c("f0", "Q", "k")
# colnames(fit_Q100_mle) <- c("f0", "Q", "k")
# # Q = 500
# fit_Q500_lp <- fit_Q500_lp %*% diag(c(1/f0, 1/Q500, 1/k))
# fit_Q500_nls <- fit_Q500_nls %*% diag(c(1/f0, 1/Q500, 1/k))
# fit_Q500_mle <- fit_Q500_mle %*% diag(c(1/f0, 1/Q500, 1/k))
# colnames(fit_Q500_lp) <- c("f0", "Q", "k")
# colnames(fit_Q500_nls) <- c("f0", "Q", "k")
# colnames(fit_Q500_mle) <- c("f0", "Q", "k")
# # merge datasets
# # Q = 1
# fit_Q1_lp <- as_tibble(fit_Q1_lp) %>% add_column(method = "LP")
# fit_Q1_nls <- as_tibble(fit_Q1_nls) %>% add_column(method = "NLS")
# fit_Q1_mle <- as_tibble(fit_Q1_mle) %>% add_column(method = "MLE")
# fit_Q1 <- bind_rows(fit_Q1_lp, fit_Q1_nls, fit_Q1_mle)
# fit_Q1 <- fit_Q1 %>% add_column(level = "Q = 1")
# # Q = 10
# fit_Q10_lp <- as_tibble(fit_Q10_lp) %>% add_column(method = "LP")
# fit_Q10_nls <- as_tibble(fit_Q10_nls) %>% add_column(method = "NLS")
# fit_Q10_mle <- as_tibble(fit_Q10_mle) %>% add_column(method = "MLE")
# fit_Q10 <- bind_rows(fit_Q10_lp, fit_Q10_nls, fit_Q10_mle) %>% 
#   add_column(level = "Q = 10")
# # Q = 100
# fit_Q100_lp <- as_tibble(fit_Q100_lp) %>% add_column(method = "LP")
# fit_Q100_nls <- as_tibble(fit_Q100_nls) %>% add_column(method = "NLS")
# fit_Q100_mle <- as_tibble(fit_Q100_mle) %>% add_column(method = "MLE")
# fit_Q100 <- bind_rows(fit_Q100_lp, fit_Q100_nls, fit_Q100_mle) %>% 
#   add_column(level = "Q = 100")
# # Q = 500
# fit_Q500_lp <- as_tibble(fit_Q500_lp) %>% add_column(method = "LP")
# fit_Q500_nls <- as_tibble(fit_Q500_nls) %>% add_column(method = "NLS")
# fit_Q500_mle <- as_tibble(fit_Q500_mle) %>% add_column(method = "MLE")
# fit_Q500 <- bind_rows(fit_Q500_lp, fit_Q500_nls, fit_Q500_mle) %>% 
#   add_column(level = "Q = 500")
# # combine all the datasets together
# fit <- bind_rows(fit_Q1, fit_Q10, fit_Q100, fit_Q500)
# # boxplot
# # Q_hat / Q
# ggplot(fit, aes(x = level, y = Q, fill = method)) + geom_boxplot()
# # k_hat / k
# ggplot(fit, aes(x = level, y = k, fill = method)) + geom_boxplot()
# # f0_hat / f0
# ggplot(fit, aes(x = level, y = f0, fill = method)) + geom_boxplot()


# nbins <- ceiling(length(f) / binSize)
# nleft <- nbins*binSize - length(f) # number of positions left
# ftmp <- c(f, rep(0, nleft)) # fill the remaining positions with 0s
# # binning the frequency domain and periodogram values
# fmat <- matrix(ftmp, nrow = binSize) # each column is a bin whose values have consecutive memory locations
# fbar <- colMeans(fmat)
# fbar[nbins] <- sum(f[(binSize*(nbins-1)+1):length(f)])/(length(f) - binSize*(nbins-1))


