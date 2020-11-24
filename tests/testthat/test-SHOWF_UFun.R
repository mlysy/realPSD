source("realPSD-testfunctions.R")

context("SHOWF UFun Tests")

test_that("The SHOWF UFun returned is the same in R and TMB", {
  ntest <- 50
  nphi <- sample(2:5,1)
  for(ii in 1:ntest) {
    # pick model
    model <- sample(c("SHOWF_log", "SHOWF_nat"), size = 1)
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    Y <- sim_Y(N)
    fs <- sim_fs()
    phi0 <- sim_showf_phi(model) # phi = c(f0, Q, Rw, Rf, alpha)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model, 
                                      method = "UFun", 
                                      f = matrix(f),
                                      # Y = matrix(Y),
                                      fs = fs),
                           parameters = list(phi = matrix(phi0)),
                           # map = map,
                           silent = TRUE, DLL = "realPSD_TMBExports")
    ufun_tmb <- function(phi) c(tmod$simulate(phi)$U)
    # check they are equal
    Phi <- replicate(nphi, sim_showf_phi(model = model))
    U_r <- apply(Phi, 2, ufun_r, f = f)
    U_tmb <- apply(Phi, 2, ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})

test_that("The SHOWF UFun (with map) returned is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5,1)
  for(ii in 1:ntest) {
    # pick model
    model <- sample(c("SHOWF_log", "SHOWF_nat"), size = 1)
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    Y <- sim_Y(N)
    fs <- sim_fs()
    phi0 <- sim_showf_phi(model) # phi = c(f0, Q, Rw, Rf, alpha)
    # create TMB model and functions
    map <- list(as.factor(c(1,2,NA,4,5)))
    phi0[3] <- 0
    tmod <- TMB::MakeADFun(data = list(model = model, 
                                      method = "UFun", 
                                      f = matrix(f),
                                      # Y = matrix(Y),
                                      fs = fs),
                           parameters = list(phi = matrix(phi0)),
                           map = map,
                           silent = TRUE, DLL = "realPSD_TMBExports")
    ufun_tmb <- function(phi) c(tmod$simulate(phi)$U)
    # check they are equal
    Phi <- replicate(nphi, sim_showf_phi(model = model))
    Phi["Rw", ] <- rep(0, nphi)
    U_r <- apply(Phi, 2, ufun_r, f = f)
    U_tmb <- apply(Phi, 2, ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})
