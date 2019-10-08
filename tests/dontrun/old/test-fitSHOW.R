# test fitSHOW
require(testthat)
require(realPSD)

source("test-functions.R")
source("fitSHOW.R")

context("fitSHOW")

test_that("Estimated parameters given by fitSHOW (MLE method) are correct", {
  ntest <- 20
  for(ii in 1:ntest) {
    # ---------- setup ----------
    # Boltzmann's constant
    Kb <- 1.381e-23
    # simualte data
    N <- sample(10:50, 1) # number of observations
    fseq <- sort(sim_f(N)) # frequency domain/sequence
    sim_exp <- rexp(N, rate = 1) # random exponential variables
    phi <- sim_phi() # phi = c(f0, f0*Q, Rw)
    # phi <- c(f0 = 33553, gamma = 3355300, Rw = 19000/1)
    fs <- sample(10:1000, 1) # sampling frequency
    k <- runif(1) # stiffness
    Temp <- sample(100:500, 1) # temperature
    bin_size <- sample(1:3, 1) # bin size
    # derive some constant inputs
    f0 <- phi[1]
    Q <- phi[2] / phi[1]
    tau <- Kb * Temp / (k * pi * phi[2])
    Aw <- phi[3] * tau
    # ---------- fitSHOW -----------
    param_fit <- fitSHOW(fseq, sim_exp, fs, f0, Q, k, Temp, Aw, bin_size,
      method = "mle")
    # ---------- TMB and R ----------
    psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion = TRUE) + Aw
    Y <- sim_exp * fs * psd
    fbar <- binning(fseq, bin_size = bin_size)
    Ybar <- binning(Y, bin_size = bin_size)
    Zbar <- log(Ybar)
    # model fitting
    obj <- TMB::MakeADFun(data = list(model = "SHOWFit",
                                       method = "MLE_nlp",
                                       f = matrix(fseq),
                                       Y = matrix(Y)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    # get tau
    get_tau <- TMB::MakeADFun(data = list(model = "SHOWFit",
                                       method = "MLE_tau",
                                       f = matrix(fseq),
                                       Y = matrix(Y)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    opt <- optim(phi, fn = obj$fn, gr = obj$gr,
               control = list(maxit = 1000))
    phi_hat <- opt$par
    tau_hat <- get_tau$fn(phi_hat)
    param_r <- rep(NA, 4)
    param_r[1] <- phi_hat[1] # f0_hat
    param_r[2] <- phi_hat[2]/phi_hat[1] # Q_hat
    param_r[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
    param_r[4] <- phi_hat[3] * tau_hat # Aw_hat
    names(param_r) <- c("f0", "Q", "k", "Aw")
    # ---------- compare (should be equal) ----------
    expect_equal(param_fit, param_r)
  }
})

test_that("Estimated parameters given by fitSHOW (LP method) are correct", {
  ntest <- 20
  for(ii in 1:ntest) {
    # ---------- setup ----------
    # Boltzmann's constant
    Kb <- 1.381e-23
    # simualte data
    N <- sample(10:50, 1) # number of observations
    fseq <- sort(sim_f(N)) # frequency domain/sequence
    sim_exp <- rexp(N, rate = 1) # random exponential variables
    phi <- sim_phi() # phi = c(f0, f0*Q, Rw)
    # phi <- c(f0 = 33553, gamma = 3355300, Rw = 19000/1)
    fs <- sample(10:1000, 1) # sampling frequency
    k <- runif(1) # stiffness
    Temp <- sample(100:500, 1) # temperature
    bin_size <- sample(1:3, 1) # bin size
    # derive some constant inputs
    f0 <- phi[1]
    Q <- phi[2] / phi[1]
    tau <- Kb * Temp / (k * pi * phi[2])
    Aw <- phi[3] * tau
    # ---------- fitSHOW -----------
    param_fit <- fitSHOW(fseq, sim_exp, fs, f0, Q, k, Temp, Aw, bin_size,
      method = "lp")
    # ---------- TMB and R ----------
    psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion = TRUE) + Aw
    Y <- sim_exp * fs * psd
    fbar <- binning(fseq, bin_size = bin_size)
    Ybar <- binning(Y, bin_size = bin_size)
    Zbar <- log(Ybar)
    # model fitting
    obj <- TMB::MakeADFun(data = list(model = "SHOWFit",
                                       method = "LP_nlp",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    # get tau
    get_tau <- TMB::MakeADFun(data = list(model = "SHOWFit",
                                       method = "LP_zeta",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    opt <- optim(phi, fn = obj$fn, gr = obj$gr,
               control = list(maxit = 1000))
    phi_hat <- opt$par
    tau_hat <- exp(get_tau$fn(phi_hat))
    param_r <- rep(NA, 4)
    param_r[1] <- phi_hat[1] # f0_hat
    param_r[2] <- phi_hat[2]/phi_hat[1] # Q_hat
    param_r[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
    param_r[4] <- phi_hat[3] * tau_hat # Aw_hat
    names(param_r) <- c("f0", "Q", "k", "Aw")
    # ---------- compare (should be equal) ----------
    expect_equal(param_fit, param_r)
  }
})

test_that("Estimated parameters given by fitSHOW (NLS method) are correct", {
  ntest <- 20
  for(ii in 1:ntest) {
    # ---------- setup ----------
    # Boltzmann's constant
    Kb <- 1.381e-23
    # simualte data
    N <- sample(10:50, 1) # number of observations
    fseq <- sort(sim_f(N)) # frequency domain/sequence
    sim_exp <- rexp(N, rate = 1) # random exponential variables
    phi <- sim_phi() # phi = c(f0, f0*Q, Rw)
    # phi <- c(f0 = 33553, gamma = 3355300, Rw = 19000/1)
    fs <- sample(10:1000, 1) # sampling frequency
    k <- runif(1) # stiffness
    Temp <- sample(100:500, 1) # temperature
    bin_size <- sample(1:3, 1) # bin size
    # derive some constant inputs
    f0 <- phi[1]
    Q <- phi[2] / phi[1]
    tau <- Kb * Temp / (k * pi * phi[2])
    Aw <- phi[3] * tau
    # ---------- fitSHOW -----------
    param_fit <- fitSHOW(fseq, sim_exp, fs, f0, Q, k, Temp, Aw, bin_size,
      method = "nls")
    # ---------- TMB and R ----------
    psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion = TRUE) + Aw
    Y <- sim_exp * fs * psd
    fbar <- binning(fseq, bin_size = bin_size)
    Ybar <- binning(Y, bin_size = bin_size)
    Zbar <- log(Ybar)
    # model fitting
    obj <- TMB::MakeADFun(data = list(model = "SHOWFit",
                                       method = "NLS_nlp",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    # get tau
    get_tau <- TMB::MakeADFun(data = list(model = "SHOWFit",
                                       method = "NLS_tau",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    opt <- optim(phi, fn = obj$fn, gr = obj$gr,
               control = list(maxit = 1000))
    phi_hat <- opt$par
    tau_hat <- get_tau$fn(phi_hat)
    param_r <- rep(NA, 4)
    param_r[1] <- phi_hat[1] # f0_hat
    param_r[2] <- phi_hat[2]/phi_hat[1] # Q_hat
    param_r[3] <- Kb * Temp / (tau_hat * pi * phi_hat[2]) # k_hat
    param_r[4] <- phi_hat[3] * tau_hat # Aw_hat
    names(param_r) <- c("f0", "Q", "k", "Aw")
    # ---------- compare (should be equal) ----------
    expect_equal(param_fit, param_r)
  }
})
