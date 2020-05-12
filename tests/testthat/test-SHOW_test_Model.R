source("realPSD-testfunctions.R")

context("SHOW-test")

test_that("The UFun returned is the same in R and TMB", {
  ntest <- 20
  nphi <- sample(2:5, 1)
  for(ii in 1:ntest) {
    # pick model
    model <- "SHOW_test"
    # ufun_r <- function(f, phi, mult_factor) {
    #   u <- phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/(phi[1]*phi[2]))^2)
    #   u <- mult_factor * u
    # }
<<<<<<< HEAD
    ufun_r <- get_ufun("SHOW_nat") # since the test model is based on SHOW_nat
=======
    ufun_r <- get_ufun("SHOW_log") # since the test model is based on SHOW_log
>>>>>>> ca74bce3fa39be73e569f1c7bd492a06bf944166
    # simulate data
    N <- sample(10:20,1)
    f <- sim_f(N)
    mult_factor <- 1 # if mult_factor == 1, then the result should be the same as test-UFunMethods.R
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model = model,
<<<<<<< HEAD
                                       method = "UFun",
                                       f = matrix(f),
                                       mult_factor = mult_factor),
=======
                                       method = "LP_nlp",
                                       f = matrix(f),
                                       mult_factor = matrix(mult_factor)),
>>>>>>> ca74bce3fa39be73e569f1c7bd492a06bf944166
                           parameters = list(phi = matrix(rep(0, 3))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
    ufun_tmb <- function(phi) c(tmod$simulate(phi)$U)
    # check they are equal
<<<<<<< HEAD
    Phi <- replicate(nphi, sim_phi(model = "SHOW_nat"))
=======
    Phi <- replicate(nphi, sim_phi(model = "SHOW_log"))
>>>>>>> ca74bce3fa39be73e569f1c7bd492a06bf944166
    U_r <- apply(Phi, 2, ufun_r, f = f)
    U_tmb <- apply(Phi, 2, ufun_tmb)
    expect_equal(U_r, U_tmb)
  }
})
