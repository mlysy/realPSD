source("realPSD-testfunctions.R")

context("MLEMethods")

test_that("MLE_ufun is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "mle"
  for(ii in 1:ntest) {
    bin_type <- "mean"
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ypsd)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$ufun()
    } else {
      scale <- mean(Ypsd) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "MLE_ufun",
                                         f = matrix(freq),
                                         Y = matrix(Ypsd/scale),
                                         scale = log(scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             ADreport = TRUE,
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    mle_ufun_tmb <- function(phi) setNames(tmod$fn(phi), NULL)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    U_r <- apply(Phi, 2, ufun_r, f = freq) # * fs
    U_tmb <- apply(Phi, 2, mle_ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})

test_that("MLE_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "mle"
  for(ii in 1:ntest) {
    bin_type <- "mean"
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ypsd)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$zeta()
    } else {
      scale <- mean(Ypsd) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "MLE_zeta",
                                         f = matrix(freq),
                                         Y = matrix(Ypsd/scale),
                                         scale = log(scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    mle_zeta_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, mle_tau_r, f = freq, Y = Ypsd,
                   ufun = ufun_r, fs = 1)
    zeta_r <- log(tau_r)
    zeta_tmb <- apply(Phi, 2, mle_zeta_tmb)
    expect_equal(zeta_r, zeta_tmb)
  }
})

test_that("MLE_nlp_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "mle"
  for(ii in 1:ntest) {
    bin_type <- "mean"
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ypsd)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nlp()
    } else {
      scale <- mean(Ypsd) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "MLE_nlp",
                                         f = matrix(freq),
                                         Y = matrix(Ypsd/scale),
                                         scale = log(scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    mle_zeta_tmb <- function(phi) tmod$simulate(phi)$zeta
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, mle_tau_r, f = freq, Y = Ypsd,
                   ufun = ufun_r, fs = 1)
    zeta_r <- log(tau_r)
    zeta_tmb <- apply(Phi, 2, mle_zeta_tmb)
    expect_equal(zeta_r, zeta_tmb)
  }
})


test_that("MLE_nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "mle"
  for(ii in 1:ntest) {
    bin_type <- "mean"
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ypsd)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nll()
    } else {
      scale <- mean(Ypsd) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "MLE_nll",
                                         f = matrix(freq),
                                         Y = matrix(Ypsd/scale),
                                         scale = log(scale)),
                             parameters = list(phi = matrix(rep(0, 3)),
                                               zeta = 0),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    mle_nll_tmb <- function(phi, zeta) tmod$fn(c(phi, zeta))
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    ## zeta <- replicate(nphi, sim_zeta())
    nll_r <- sapply(1:nphi, function(ii) {
      mle_nll_r(phi = Phi[,ii], tau = exp(zeta[ii]), Y = Ypsd,
        f = freq, ufun = ufun_r, fs = 1)
    })
    nll_tmb <- sapply(1:nphi, function(ii) mle_nll_tmb(Phi[,ii], zeta[ii]))
    expect_equal(nll_r, nll_tmb)
  }
})

test_that("MLE_nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "mle"
  for(ii in 1:ntest) {
    bin_type <- "mean"
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ypsd)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nlp()
    } else {
      scale <- mean(Ypsd) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "MLE_nlp",
                                         f = matrix(freq),
                                         Y = matrix(Ypsd/scale),
                                         scale = log(scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    mle_nlp_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      mle_nlp_r(phi = phi, Y = Ypsd, f = freq, ufun = ufun_r, fs = 1)
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
