## tmb_recompile()

## require(realPSD)
## require(TMB)
## require(testthat)
source("realPSD-testfunctions.R")

context("LPMethods")

test_that("LP_ufun is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "lp"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    # create TMB model and functions
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Zbar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$ufun()
    } else {
      scale <- mean(Zbar) * runif(1, .9, 1.1)
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "LP_ufun",
                                         fbar = matrix(fbar),
                                         Zbar = matrix(Zbar - scale),
                                         scale = c(scale - bin_loc, bin_scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             ADreport = TRUE,
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    lp_ufun_tmb <- function(phi) setNames(tmod$fn(phi), NULL)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    U_r <- apply(Phi, 2, ufun_r, f = fbar) # * fs
    U_tmb <- apply(Phi, 2, lp_ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})

test_that("LP_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "lp"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Zbar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                est_type = est_type,
                                bin_size = B, bin_type = bin_type)
      tmod <- test_obj$zeta()
    } else {
      scale <- mean(Zbar) * runif(1, .9, 1.1)
      # create TMB model and functions
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "LP_zeta",
                                         fbar = matrix(fbar),
                                         Zbar = matrix(Zbar - scale),
                                         scale = c(scale - bin_loc, bin_scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    lp_zeta_tmb <- function(phi) tmod$fn(phi) # objective function value returned by TMB
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    zeta_r <- apply(Phi, 2, lp_zeta_r, fbar = fbar, Zbar = Zbar,
                    ufun = ufun_r, B = B, fs = 1, bin_type = bin_type)
    zeta_tmb <- apply(Phi, 2, lp_zeta_tmb)
    expect_equal(zeta_r, zeta_tmb)
  }
})

test_that("LP_nlp_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "lp"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Zbar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nlp()
    } else {
      scale <- mean(Zbar) * runif(1, .9, 1.1)
      # create TMB model and functions
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "LP_nlp",
                                         fbar = matrix(fbar),
                                         Zbar = matrix(Zbar - scale),
                                         scale = c(scale - bin_loc, bin_scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    lp_zeta_tmb <- function(phi) tmod$simulate(phi)$zeta
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    zeta_r <- apply(Phi, 2, lp_zeta_r, fbar = fbar, Zbar = Zbar,
                    ufun = ufun_r, fs = 1, B = B, bin_type = bin_type)
    zeta_tmb <- apply(Phi, 2, lp_zeta_tmb)
    expect_equal(zeta_r, zeta_tmb)
  }
})

test_that("LP_nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "lp"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Zbar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nll()
    } else {
      scale <- mean(Zbar) * runif(1, .9, 1.1)
      # create TMB model and functions
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "LP_nll",
                                         fbar = matrix(fbar),
                                         Zbar = matrix(Zbar - scale),
                                         scale = c(scale - bin_loc, bin_scale)),
                             parameters = list(phi = matrix(rep(0, 3)),
                                               zeta = 0),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    lp_nll_tmb <- function(phi, zeta) tmod$fn(c(phi, zeta))
    # check they are equal
    nll_r <- sapply(1:nphi, function(ii) {
      lp_nll_r(phi = Phi[,ii], zeta = zeta[ii],
               Zbar = Zbar, fbar = fbar, ufun = ufun_r, fs = 1, B = B, bin_type = bin_type)
    })
    nll_tmb <- sapply(1:nphi, function(ii) lp_nll_tmb(Phi[,ii], zeta[ii]))
    expect_equal(nll_r, nll_tmb)
  }
})

test_that("LP_nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "lp"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Zbar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$nlp()
    } else {
      scale <- mean(Zbar) * runif(1, .9, 1.1)
      # create TMB model and functions
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "LP_nlp",
                                         fbar = matrix(fbar),
                                         Zbar = matrix(Zbar - scale),
                                         scale = c(scale - bin_loc, bin_scale)),
                             parameters = list(phi = matrix(rep(0, 3))),
                             silent = TRUE, DLL = "realPSD_TMBExports")
    }
    lp_nlp_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      lp_nlp_r(phi = phi, Zbar = Zbar, fbar = fbar,
               ufun = ufun_r, fs = 1, B = B, bin_type = bin_type)
    })
    nlp_tmb <- apply(Phi, 2, lp_nlp_tmb)
    expect_equal(nlp_r, nlp_tmb)
  }
})

test_that("LP_res is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  est_type <- "lp"
  for(ii in 1:ntest) {
    bin_type <- sample(c("mean", "median"), 1)
    with_r6 <- sample(c(TRUE, FALSE), 1)
    list2env(sim_setup(nphi, est_type = est_type, bin_type = bin_type),
             envir = environment())
    if(with_r6) {
      scale <- mean(Zbar)
      test_obj <- test_model$new(freq = freq, Ypsd = Ypsd,
                                 est_type = est_type,
                                 bin_size = B, bin_type = bin_type)
      tmod <- test_obj$resid()
    } else {
      scale <- mean(Zbar) * runif(1, .9, 1.1)
      # create TMB model and functions
      tmod <- TMB::MakeADFun(data = list(model = model,
                                         method = "LP_res",
                                         fbar = matrix(fbar),
                                         Zbar = matrix(Zbar - scale),
                                         ## fs = fs,
                                         scale = c(scale - bin_loc, bin_scale)),
                             parameters = list(phi = matrix(rep(0, 3)), zeta = 0),
                             silent = TRUE,
                             ADreport = TRUE,
                             DLL = "realPSD_TMBExports")
    }
    lp_res_tmb <- function(phi, zeta) setNames(tmod$fn(c(phi, zeta)), nm = NULL)
    # check they are equal
    ## Phi <- replicate(nphi, sim_phi())
    ## zeta <- replicate(nphi, sim_zeta())
    res_r <- sapply(1:nphi, function(ii) {
      lp_res_r(phi = Phi[,ii], zeta = zeta[ii], fbar = fbar, Zbar = Zbar,
               ufun = ufun_r, fs = 1, B = B, bin_type = bin_type)
    })
    res_tmb <- sapply(1:nphi, function(ii) {
      lp_res_tmb(phi = Phi[,ii], zeta = zeta[ii])
    })
    expect_equal(res_r, res_tmb)
  }
})


#--- scratch -------------------------------------------------------------------

## N <- 10
## fbar <- sim_f(N)
## ## fbar[1] <- 1e-200
## ## fbar <- c(rep(0, N-1),1)
## Zbar <- sim_Zbar(N)
## phi <- c(1, 1, 0)
## zeta <- 0
## ufun_r(fbar, phi)
## nllik_r(phi, zeta, Zbar, fbar)
## sum(Zbar^2)

## tmod <- TMB::MakeADFun(data = list(model = "SHOWFit",
##                                    method = "nll",
##                                    fbar = matrix(fbar),
##                                    Zbar = matrix(Zbar)),
##                        parameters = list(phi = matrix(0, 3, 1),
##                                          zeta = 0),
##                        silent = TRUE, DLL = "realPSD_TMBExports")

## phi <- sim_phi()
## zeta <- sim_zeta()
## ## logUbar <- log(ufun_r(fbar, phi))
## ## tmod$report(c(phi, zeta))
## tmod$fn(c(phi, zeta))
## ## fbar[1]^2/phi[1]
## ## sum((Zbar - log(fbar^2/phi[1]^2 + 1) - zeta)^2)
## ## logUbar[1]
## nllik_r(phi, zeta, Zbar, fbar)



## nllik_tmb <- function(phi, zeta) tmod$fn(c(phi, zeta))


## N <- sample(10:20,1)
## fbar <- sim_f(N)
## Zbar <- sim_Zbar(N)
## phi <- sim_phi()
## zeta_hat <- zeta_r(fbar, Zbar, phi)

## zeta_seq <- seq(zeta_hat - 3 * abs(zeta_hat), zeta_hat + 3 * abs(zeta_hat),
##                 len = 100)
## zeta_ll <- sapply(zeta_seq, nllik_r, phi = phi, Zbar = Zbar, fbar = fbar)
## plot(zeta_seq, zeta_ll, type = "l")
## abline(v = zeta_hat)
