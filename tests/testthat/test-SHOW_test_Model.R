source("realPSD-testfunctions.R")

context("SHOW-test")

test_that("The UFun returned is the same in R and TMB", {
  ntest <- 20
  for(ii in 1:ntest) {
    # pick model
    model <- "SHOW_test"
    # simulate data
    N <- sample(10:20,1)
    f <- matrix(sim_f(N))
    alpha <- rexp(1) # multiplicative factor
    phi <- matrix(rexp(3))
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model, f = f),
                           parameters = list(phi = matrix(rep(0, 3)), alpha = 0),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    expect_equal(
      tmod$simulate(c(phi, alpha))$U/tmod$simulate(c(phi, 1))$U, 
      matrix(rep(alpha, N), ncol = 1)
    )
  }
})
