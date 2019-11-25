source("realPSD-testfunctions.R")

context("MLEMethods")

test_that("MLE_tau is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    Y <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "MLE_tau",
                                       f = matrix(f),
                                       Y = matrix(Y),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    mle_tau_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, mle_tau_r, f = f, Y = Y,
                    ufun = ufun_r, fs = fs)
    tau_tmb <- apply(Phi, 2, mle_tau_tmb)
    expect_equal(tau_r, tau_tmb)
  }
})

test_that("MLE_nlp_tau is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    Y <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "MLE_nlp",
                                       f = matrix(f),
                                       Y = matrix(Y),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    mle_tau_tmb <- function(phi) tmod$simulate(phi)$tau
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, mle_tau_r, f = f, Y = Y,
                    ufun = ufun_r, fs = fs)
    tau_tmb <- apply(Phi, 2, mle_tau_tmb)
    expect_equal(tau_r, tau_tmb)
  }
})


test_that("MLE_nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    Y <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "MLE_nll",
                                       f = matrix(f),
                                       Y = matrix(Y),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3)),
                                             tau = 0),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    mle_nll_tmb <- function(phi, tau) tmod$fn(c(phi, tau))
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    tau <- replicate(nphi, sim_tau())
    nll_r <- sapply(1:nphi, function(ii) {
      mle_nll_r(phi = Phi[,ii], tau = tau[ii], Y = Y,
        f = f, ufun = ufun_r, fs = fs)
    })
    nll_tmb <- sapply(1:nphi, function(ii) mle_nll_tmb(Phi[,ii], tau[ii]))
    expect_equal(nll_r, nll_tmb)
  }
})

test_that("MLE_nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    Y <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "MLE_nlp",
                                       f = matrix(f),
                                       Y = matrix(Y),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    mle_nlp_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      mle_nlp_r(phi = phi, Y = Y, f = f, ufun = ufun_r, fs = fs)
    })
    nlp_tmb <- apply(Phi, 2, mle_nlp_tmb)
    expect_equal(nlp_r, nlp_tmb)
  }
})

#--- scratch -------------------------------------------------------------------

## # SHOW normalized PSD
## ufun_r <- function(f, phi) {
##   phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/phi[2])^2)
## }

## # full negative loglikelihood
## nllik_r <- function(phi, tau, Y, f) {
##   U <- ufun_r(f, phi)
##   S <- tau * U
##   if(any(S < 0)) browser()
##   sum(Y/S + log(S))
## }

## # conditional maximum
## tau_r <- function(f, Y, phi) {
##   U <- ufun_r(f, phi)
##   mean(Y/U)
## }

## # profile negative loglikelihood
## nlprof_r <- function(phi, Y, f) {
##   U <- ufun_r(f, phi)
##   tau <- tau_r(f, Y, phi)
##   length(Y) * (1 + log(tau)) + sum(log(U))
## }

## sim_f <- function(n) runif(n, 0, 2*n)
## sim_Y <- function(n) rchisq(n, df = 2)
## sim_phi <- function() c(f0 = runif(1, 100, 1000),
##                         gamma = rexp(1), Rw = rexp(1))
## sim_tau <- function() rexp(1)

## n <- 50
## f <- sim_f(n)
## Y <- sim_Y(n)
## phi <- sim_phi()
## tau_hat <- tau_r(f, Y, phi)
## tau_seq <- seq(tau_hat - .5 * abs(tau_hat),
##                tau_hat + .5 * abs(tau_hat), len = 100)
## tau_nll <- sapply(tau_seq, function(tau) nllik_r(phi, tau, Y, f))
## ## plot(tau_seq, tau_nll, type = "l")
## ## abline(v = tau_hat)

## n <- 50
## f <- sim_f(n)
## Y <- sim_Y(n)
## phi <- sim_phi()
## tau_hat <- tau_r(f, Y, phi)
## nllik_r(phi, tau_hat, Y, f)
## nlprof_r(phi, Y, f)
