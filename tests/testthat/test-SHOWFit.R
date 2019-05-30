## tmb_recompile()

## require(realPSD)
## require(TMB)
## require(testthat)

context("SHOWFit")

source("realPSD-testfunctions.R")

test_that("UFun is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5,1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "UFun",
                                       f = matrix(f)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    ufun_tmb <- function(phi) c(tmod$simulate(phi)$U)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    U_r <- apply(Phi, 2, ufun_r, f = f)
    U_tmb <- apply(Phi, 2, ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})

test_that("zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "zeta",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    zeta_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    z_r <- apply(Phi, 2, zeta_r, fbar = fbar, Zbar = Zbar)
    z_tmb <- apply(Phi, 2, zeta_tmb)
    expect_equal(z_r, z_tmb)
  }
})

test_that("nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "nll",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(rep(0, 3)),
                                             zeta = 0),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nllik_tmb <- function(phi, zeta) tmod$fn(c(phi, zeta))
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    zeta <- replicate(nphi, sim_zeta())
    nll_r <- sapply(1:nphi, function(ii) {
      nllik_r(Phi[,ii], zeta[ii], Zbar, fbar)
    })
    nll_tmb <- sapply(1:nphi, function(ii) nllik_tmb(Phi[,ii], zeta[ii]))
    expect_equal(nll_r, nll_tmb)
  }
})

test_that("nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "nlp",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nlprof_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      zeta <- zeta_r(fbar, Zbar, phi)
      nllik_r(phi, zeta, Zbar, fbar)
    })
    nlp_tmb <- apply(Phi, 2, nlprof_tmb)
    expect_equal(nlp_r, nlp_tmb)
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

## tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
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
