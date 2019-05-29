require(realPSD)
require(TMB)
require(testthat)

context("SHOWFit")

sim_f <- function(n) runif(n, 0, 2*n)
sim_Zbar <- function(n) rnorm(n)
sim_phi <- function() c(f0 = runif(1, 100, 1000), gamma = rexp(1), Rw = rexp(1))
ufun_r <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/phi[2])^2)
}
zeta_r <- function(fbar, Zbar, phi) {
  mean(Zbar - log(ufun_r(fbar, phi)))
}
lpobj_r <- function(phi, zeta, Zbar, fbar) {
  logUbar <- log(ufun_r(fbar, phi))
  sum((Zbar - zeta - logUbar)^2)
}

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

#--- scratch -------------------------------------------------------------------

## N <- sample(10:20,1)
## fbar <- sim_f(N)
## Zbar <- sim_Zbar(N)
## phi <- sim_phi()
## zeta_hat <- zeta_r(fbar, Zbar, phi)

## zeta_seq <- seq(zeta_hat - 3 * abs(zeta_hat), zeta_hat + 3 * abs(zeta_hat),
##                 len = 100)
## zeta_ll <- sapply(zeta_seq, lpobj_r, phi = phi, Zbar = Zbar, fbar = fbar)
## plot(zeta_seq, zeta_ll, type = "l")
## abline(v = zeta_hat)
