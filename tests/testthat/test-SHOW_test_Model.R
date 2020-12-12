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
    tmod1 <- TMB::MakeADFun(data = list(model = model, method = "UFun",
                                        f = f, alpha = alpha),
                            parameters = list(phi = matrix(rep(0, 3))),
                            ADreport = TRUE,
                            silent = TRUE, DLL = "realPSD_TMBExports")
    tmod2 <- TMB::MakeADFun(data = list(model = model, method = "UFun",
                                        f = f, alpha = 1),
                            parameters = list(phi = matrix(rep(0, 3))),
                            ADreport = TRUE,
                            silent = TRUE, DLL = "realPSD_TMBExports")
    expect_equal(
      setNames(tmod1$fn(phi)/tmod2$fn(phi), NULL),
      rep(alpha, N)
    )
  }
})
