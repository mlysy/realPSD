#--- initialization of fractional OU model object ---------------------------------------------

require(TMB)
require(testthat)

model <- "fou"
compile(paste0(model, ".cpp"))
dyn.load(dynlib(model))

# helper function: calculate the psd in R
fou_psd_r <- function(f, phi) {
 psd <- abs(f)^(1-2*phi[1]) / (f^2 + phi[2]^2)
}

context("Test: fractional OU model")

test_that("The UFun returned by TMB (C++) is the same as that given by R", {
  ntests <- 20
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    f <- matrix(runif(N, 0, 2*N))
    phi <- matrix(c(runif(1), rexp(1))) # phi = c(H, gamma) where 0 < H < 1, gamma > 0
    obj <- MakeADFun(data = list(f = f),
                 parameters = list(phi = matrix(c(runif(1), rexp(1)))),
                 DLL = model, silent = TRUE)
    psd_tmb <- obj$simulate(c(phi))$U
    psd_r <- fou_psd_r(f, phi)
    expect_equal(psd_tmb, psd_r)
  }
})

test_that("sigma^2 returned by TMB (C++) is the same as that given by R", {
  ntests <- 20
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    f <- matrix(runif(N, 0, 2*N))
    phi <- matrix(c(runif(1), rexp(1))) 
    obj <- MakeADFun(data = list(f = f), 
      parameters = list(phi = matrix(0.5, 2)), 
      DLL = model, silent = TRUE)
  }
})

