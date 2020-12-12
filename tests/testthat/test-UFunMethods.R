context("UFunMethods")

source("realPSD-testfunctions.R")

test_that("UFun is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5,1)
  for(ii in 1:ntest) {
    # pick model
    model <- sim_model()
    ufun_r <- get_ufun(model)
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "UFun",
                                       f = matrix(f)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           ADreport = TRUE,
                           silent = TRUE, DLL = "realPSD_TMBExports")
    ufun_tmb <- function(phi) setNames(tmod$fn(phi), NULL)
    # check they are equal
    Phi <- replicate(nphi, sim_phi(model = model))
    U_r <- apply(Phi, 2, ufun_r, f = f)
    U_tmb <- apply(Phi, 2, ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})
