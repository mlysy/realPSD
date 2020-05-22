source("realPSD-testfunctions.R")

context("SHOW-test")

test_that("The UFun returned is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- "SHOW_test"
    ufun_r <- get_ufun("SHOW_log") # since the test model is based on SHOW_log
    # simulate data
    N <- sample(10:20,1)
    fbar <- sim_f(N)
    Zbar <- sim_Zbar(N)
    fs <- sim_fs()
    mult_factor <- 1 # if mult_factor == 1, then the result should be the same as test-UFunMethods.R
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
                                       method = "LP_nlp",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar),
                                       fs = fs,
                                       mult_factor = matrix(mult_factor)),
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    ufun_tmb <- function(phi) c(tmod$simulate(phi)$U)
    # check they are equal
    Phi <- replicate(nphi, sim_phi(model = "SHOW_log"))
    U_r <- apply(Phi, 2, ufun_r, f = f)
    U_tmb <- apply(Phi, 2, ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})
