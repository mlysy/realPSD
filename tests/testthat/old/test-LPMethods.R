## tmb_recompile()

## require(realPSD)
## require(TMB)
## require(testthat)

context("LPMethods")

source("realPSD-testfunctions.R")

test_that("LP_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "LP_zeta",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    lp_zeta_tmb <- function(phi) tmod$fn(phi) # objective function value returned by TMB
    lp_zeta_gr_tmb <- function(phi) tmod$gr(phi)[1,] # gradient returned by TMB
    lp_zeta_he_tmb <- function(phi) tmod$he(phi) # Hessian returned by TMB
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    # fn
    zeta_r <- apply(Phi, 2, lp_zeta_r, fbar = fbar, Zbar = Zbar,
                    ufun = ufun_r, fs = fs)
    zeta_tmb <- apply(Phi, 2, lp_zeta_tmb)
    # gr
    zeta_gr_r <- apply(Phi, 2, lp_zeta_gr_r, fbar = fbar, Zbar = Zbar,
                    ufun = ufun_r, fs = fs)
    zeta_gr_tmb <- apply(Phi, 2, lp_zeta_gr_tmb)
    # hessian
    zeta_he_r <- apply(Phi, 2, lp_zeta_he_r, fbar = fbar, Zbar = Zbar,
<<<<<<< HEAD:tests/testthat/test-LPMethods.R
                    ufun = ufun_r, fs = fs) 
    zeta_he_tmb <- apply(Phi, 2, lp_zeta_he_tmb)
    # final check
    expect_equal(zeta_r, zeta_tmb)
    expect_equal(zeta_gr_r, zeta_gr_tmb, tolerance=1e-5)
    expect_equal(zeta_he_r, zeta_he_tmb, tolerance=1e-5)
=======
                    ufun = ufun_r, fs = fs)
    zeta_he_tmb <- apply(Phi, 2, lp_zeta_he_tmb)
    # final check
    expect_equal(zeta_r, zeta_tmb)
    ## expect_equal(zeta_gr_r, zeta_gr_tmb, tolerance=1e-5)
    ## expect_equal(zeta_he_r, zeta_he_tmb, tolerance=1e-5)
>>>>>>> ca74bce3fa39be73e569f1c7bd492a06bf944166:tests/testthat/old/test-LPMethods.R
  }
})

test_that("LP_nlp_zeta is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "LP_nlp",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    lp_zeta_tmb <- function(phi) tmod$simulate(phi)$zeta
    # lp_zeta_gr_tmb <- function(phi) tmod$gr(phi)[1,] # gradient returned by TMB
    # lp_zeta_he_tmb <- function(phi) tmod$he(phi) # Hessian returned by TMB
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    # fn
    zeta_r <- apply(Phi, 2, lp_zeta_r, fbar = fbar, Zbar = Zbar,
                    ufun = ufun_r, fs = fs)
    zeta_tmb <- apply(Phi, 2, lp_zeta_tmb)
    # # gr
    # zeta_gr_r <- apply(Phi, 2, lp_zeta_gr_r, fbar = fbar, Zbar = Zbar,
    #                 ufun = ufun_r, fs = fs)
    # zeta_gr_tmb <- apply(Phi, 2, lp_zeta_gr_tmb)
    # # hessian
    # zeta_he_r <- apply(Phi, 2, lp_zeta_he_r, fbar = fbar, Zbar = Zbar,
    #                 ufun = ufun_r, fs = fs)
    # zeta_he_tmb <- apply(Phi, 2, lp_zeta_he_tmb)
    expect_equal(zeta_r, zeta_tmb)
    # expect_equal(zeta_gr_r, zeta_gr_tmb, tolerance=1e-5)
    # expect_equal(zeta_he_r, zeta_he_tmb, tolerance=1e-5)
  }
})

test_that("LP_nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "LP_nll",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3)),
                                             zeta = 0),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    lp_nll_tmb <- function(phi, zeta) tmod$fn(c(phi, zeta))
    lp_nll_gr_tmb <- function(phi) tmod$gr(c(phi, zeta))[1,] # gradient returned by TMB
    lp_nll_he_tmb <- function(phi) tmod$he(c(phi, zeta)) # Hessian returned by TMB
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    zeta <- replicate(nphi, sim_zeta())
    nll_r <- sapply(1:nphi, function(ii) {
      lp_nll_r(phi = Phi[,ii], zeta = zeta[ii],
               Zbar = Zbar, fbar = fbar, ufun = ufun_r, fs = fs)
    })
    nll_tmb <- sapply(1:nphi, function(ii) lp_nll_tmb(Phi[,ii], zeta[ii]))
    # gr
    nll_gr_r <- sapply(1:nphi, function(ii) {
      lp_nll_gr_r(phi = Phi[,ii], zeta = zeta[ii],
               Zbar = Zbar, fbar = fbar, ufun = ufun_r, fs = fs)
    })
    ## nll_gr_tmb <- sapply(1:nphi, function(ii) lp_nll_gr_tmb(Phi[,ii], zeta[ii]))
    # hessian
    nll_he_r <- sapply(1:nphi, function(ii) {
      lp_nll_he_r(phi = Phi[,ii], zeta = zeta[ii],
               Zbar = Zbar, fbar = fbar, ufun = ufun_r, fs = fs)
    })
    ## nll_he_tmb <- sapply(1:nphi, function(ii) lp_nll_he_tmb(Phi[,ii], zeta[ii]))
    expect_equal(nll_r, nll_tmb)
    ## expect_equal(nll_gr_r, nll_gr_tmb, tolerance=1e-5)
    ## expect_equal(nll_he_r, nll_he_tmb, tolerance=1e-5)
  }
})

test_that("LP_nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "LP_nlp",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    lp_nlp_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      lp_nlp_r(phi = phi, Zbar = Zbar, fbar = fbar, ufun = ufun_r, fs = fs)
      ## zeta <- zeta_r(fbar, Zbar, phi)
      ## nllik_r(phi, zeta, Zbar, fbar)
    })
    nlp_tmb <- apply(Phi, 2, lp_nlp_tmb)
    expect_equal(nlp_r, nlp_tmb)
  }
})

test_that("LP_res is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "LP_res",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE,
                           ADreport = TRUE,
                           DLL = "realPSD_TMBExports")
    lp_res_tmb <- function(phi) setNames(tmod$fn(phi), nm = NULL)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    res_r <- apply(Phi, 2, lp_res_r, fbar = fbar, Zbar = Zbar,
                    ufun = ufun_r, fs = fs)
    res_tmb <- apply(Phi, 2, lp_res_tmb)
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
