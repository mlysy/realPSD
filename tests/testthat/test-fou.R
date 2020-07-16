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
    obj <- MakeADFun(data = list(f = f, method = "UFun"),
                 parameters = list(phi = matrix(c(runif(1), rexp(1)))),
                 DLL = model, silent = TRUE)
    psd_tmb <- obj$simulate(c(phi))$U
    psd_r <- fou_psd_r(f, phi)
    expect_equal(psd_tmb, psd_r)
  }
})

# helpfer function
mle_tau_r <- function(f, Y, phi, ufun, fs) {
  U <- fs * ufun(f, phi)
  mean(Y/U)
}
sim_Y <- function(n) rchisq(n, df = 2)
# normalized PSD
# phi = c(H, gamma)
fou_ufun <- function(f, phi) {
  abs(f)^(1-2*phi[1]) / (f^2 + phi[2]^2)
}
# simulate the sampling freq
sim_fs <- function() sample(10:1000, size = 1)

test_that("sigma^2 returned by TMB (C++) is the same as that given by R", {
  ntests <- 20
  for(ii in 1:ntests) {
    N <- sample(50:100, 1)
    f <- matrix(runif(N, 0, 2*N))
    fs <- sim_fs()
    phi <- matrix(c(runif(1), rexp(1))) 
    Y <- matrix(sim_Y(N))
    # TMB obj
    obj <- MakeADFun(
      data = list(f = f, Y = Y, fs = fs, method = "MLE_tau"), 
      parameters = list(phi = matrix(0.5, 2)), 
      DLL = model, silent = TRUE
    )
    # define mle_tau from TMB
    mle_tau_tmb <- function(phi) obj$fn(phi)
    # comparison
    tau_tmb <- mle_tau_tmb(phi)
    tau_r <- mle_tau_r(f, Y, phi, ufun = fou_ufun, fs)
    expect_equal(tau_tmb, tau_r)
  }
})

