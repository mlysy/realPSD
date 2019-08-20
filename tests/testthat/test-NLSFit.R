source("realPSD-testfunctions.R")

context("NLSFit")

test_that("NLS_tau is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Ybar <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "NLS_tau",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nls_tau_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    tau_r <- apply(Phi, 2, nls_tau_r, fbar = fbar, Ybar = Ybar,
                    ufun = show_ufun, fs = fs)
    tau_tmb <- apply(Phi, 2, nls_tau_tmb)
    expect_equal(tau_r, tau_tmb)
  }
})

test_that("NLS_nll is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Ybar <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "NLS_nll",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3)),
                                             tau = 0),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nls_nll_tmb <- function(phi, tau) tmod$fn(c(phi, tau))
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    tau <- replicate(nphi, sim_tau())
    nll_r <- sapply(1:nphi, function(ii) {
      nls_nll_r(phi = Phi[,ii], tau = tau[ii], Ybar = Ybar, 
        fbar = fbar, ufun = show_ufun, fs = fs)
    })
    nll_tmb <- sapply(1:nphi, function(ii) nls_nll_tmb(Phi[,ii], tau[ii]))
    expect_equal(nll_r, nll_tmb)
  }
})

test_that("NLS_res is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Ybar <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "NLS_res",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3)),
                                             tau = 0),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nls_res_tmb <- function(phi, tau) c(tmod$simulate(c(phi, tau))$RES)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    tau <- replicate(nphi, sim_tau())
    res_r <- sapply(1:nphi, function(ii) {
      nls_res_r(phi = Phi[,ii], tau = tau[ii], Ybar = Ybar, 
        fbar = fbar, ufun = show_ufun, fs = fs)
    })
    res_tmb <- sapply(1:nphi, function(ii) nls_res_tmb(Phi[,ii], tau[ii]))
    expect_equal(res_r, res_tmb)
  }
})

test_that("NLS_nlp is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Ybar <- sim_Y(N)
    fs <- sim_fs()
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "NLS_nlp",
                                       fbar = matrix(fbar),
                                       Ybar = matrix(Ybar),
                                       fs = fs),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    nls_nlp_tmb <- function(phi) tmod$fn(phi)
    # check they are equal
    Phi <- replicate(nphi, sim_phi())
    nlp_r <- apply(Phi, 2, function(phi) {
      nls_nlp_r(phi = phi, Ybar = Ybar, fbar = fbar, ufun = show_ufun, fs = fs)
    })
    nlp_tmb <- apply(Phi, 2, nls_nlp_tmb)
    expect_equal(nlp_r, nlp_tmb)
  }
})

#--- scratch -------------------------------------------------------------------

## n <- 10
## fbar <- sim_f(n)
## Ybar <- sim_Y(n)
## phi <- sim_phi()
## tau_hat <- nls_tau_r(fbar = fbar, Ybar = Ybar,
##                      phi = phi, ufun = show_ufun)
## tau_seq <- seq(tau_hat - .1 * abs(tau_hat),
##                tau_hat + .1 * abs(tau_hat), len = 100)
## tau_nll <- sapply(tau_seq, function(tau) {
##   nls_nll_r(phi = phi, tau = tau,
##             Ybar = Ybar, fbar = fbar, ufun = show_ufun)
## })
## plot(tau_seq, tau_nll, type = "l")
## abline(v = tau_hat)

