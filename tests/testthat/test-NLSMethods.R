source("realPSD-testfunctions.R")

context("NLSMethods")

test_that("NLS_ufun is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "nls"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ybar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$ufun()
    } else {
      scale <- mean(Ybar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "NLS_ufun",
                                         fbar = matrix(fbar),
                                         Ybar = matrix(Ybar/scale),
                                         scale = log(scale) - log(bin_loc)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             ADreport = TRUE,
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    nls_ufun_tmb <- function(phi) setNames(tmod$fn(phi), NULL)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    U_r <- apply(Phi, 2, ufun_r, f = fbar) # * fs
    U_tmb <- apply(Phi, 2, nls_ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})

test_that("NLS_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "nls"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ybar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$zeta()
    } else {
      scale <- mean(Ybar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "NLS_zeta",
                                         fbar = matrix(fbar),
                                         Ybar = matrix(Ybar/scale),
                                         scale = log(scale) - log(bin_loc)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    nls_zeta_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, nls_tau_r, fbar = fbar, Ybar = Ybar,
                   ufun = ufun_r, fs = 1, bin_type = bin_type)
    zeta_r <- log(tau_r)
    zeta_tmb <- apply(Phi, 2, nls_zeta_tmb)
    expect_equal(zeta_r, zeta_tmb)
  }
})

test_that("NLS_nlp_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "nls"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ybar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nlp()
    } else {
      scale <- mean(Ybar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "NLS_nlp",
                                         fbar = matrix(fbar),
                                         Ybar = matrix(Ybar/scale),
                                         scale = log(scale) - log(bin_loc)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    nls_zeta_tmb <- function(phi) tmod$simulate(phi)$zeta
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, nls_tau_r, fbar = fbar, Ybar = Ybar,
                   ufun = ufun_r, fs = 1, bin_type = bin_type)
    zeta_r <- log(tau_r)
    zeta_tmb <- apply(Phi, 2, nls_zeta_tmb)
    expect_equal(zeta_r, zeta_tmb)
  }
})

test_that("NLS_nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "nls"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ybar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nll()
    } else {
      scale <- mean(Ybar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "NLS_nll",
                                         fbar = matrix(fbar),
                                         Ybar = matrix(Ybar/scale),
                                         scale = log(scale) - log(bin_loc)),
                             parameters = list(phi = matrix(rep(0, 3)),
                                               zeta = 0),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    nls_nll_tmb <- function(phi, zeta) tmod$fn(c(phi, zeta))
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    ## zeta <- replicate(nphi, sim_zeta()) + log(scale)
    nll_r <- sapply(1:nphi, function(ii) {
      nls_nll_r(phi = Phi[,ii], tau = exp(zeta[ii]), Ybar = Ybar,
        fbar = fbar, ufun = ufun_r, fs = 1, bin_type = bin_type)
    })
    nll_tmb <- sapply(1:nphi, function(ii) nls_nll_tmb(Phi[,ii], zeta[ii]))
    expect_equal(nll_r, nll_tmb * scale^2)
  }
})

test_that("NLS_nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "nls"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ybar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nlp()
    } else {
      scale <- mean(Ybar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "NLS_nlp",
                                         fbar = matrix(fbar),
                                         Ybar = matrix(Ybar/scale),
                                         scale = log(scale) - log(bin_loc)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    nls_nlp_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      nls_nlp_r(phi = phi, Ybar = Ybar, fbar = fbar, ufun = ufun_r,
                fs = 1, bin_type = bin_type)
    })
    nlp_tmb <- apply(Phi, 2, nls_nlp_tmb)
    expect_equal(nlp_r, nlp_tmb * scale^2)
  }
})

test_that("NLS_res is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "nls"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Ybar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$resid()
    } else {
      scale <- mean(Ybar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "NLS_res",
                                         fbar = matrix(fbar),
                                         Ybar = matrix(Ybar/scale),
                                         scale = log(scale) - log(bin_loc)),
                             parameters = list(phi = matrix(rep(0, 3)), zeta = 0),
                             silent = TRUE, ADreport = TRUE,
                             DLL = "realPSD_TMBExports")
    }
    nls_res_tmb <- function(phi, tau) setNames(tmod$fn(c(phi, tau)), nm = NULL)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    ## zeta <- replicate(nphi, sim_zeta()) + log(scale)
    ## res_r <- sapply(Phi, 2, nls_res_r, fbar = fbar, Ybar = Ybar,
    ##                 ufun = ufun_r, fs = fs)
    ## res_tmb <- apply(Phi, 2, nls_res_tmb)
    res_r <- sapply(1:nphi, function(ii) {
      nls_res_r(phi = Phi[,ii], tau = exp(zeta[ii]), Ybar = Ybar,
        fbar = fbar, ufun = ufun_r, fs = 1, bin_type = bin_type)
    })
    res_tmb <- sapply(1:nphi, function(ii) nls_res_tmb(Phi[,ii], zeta[ii]))
    expect_equal(res_r, res_tmb * scale)
  }
})

#--- scratch -------------------------------------------------------------------

## n <- 10
## fbar <- sim_f(n)
## Ybar <- sim_Y(n)
## phi <- sim_phi()
## tau_hat <- nls_tau_r(fbar = fbar, Ybar = Ybar,
##                      phi = phi, ufun = ufun_r)
## tau_seq <- seq(tau_hat - .1 * abs(tau_hat),
##                tau_hat + .1 * abs(tau_hat), len = 100)
## tau_nll <- sapply(tau_seq, function(tau) {
##   nls_nll_r(phi = phi, tau = tau,
##             Ybar = Ybar, fbar = fbar, ufun = ufun_r)
## })
## plot(tau_seq, tau_nll, type = "l")
## abline(v = tau_hat)

