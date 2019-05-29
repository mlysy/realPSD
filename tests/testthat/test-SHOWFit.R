## require(realPSD)
## require(TMB)

context("SHOWFit")

sim_f <- function(n) runif(n, 0, 2*n)
sim_phi <- function() c(f0 = runif(1, 100, 1000), gamma = rexp(1), Rw = rexp(1))
ufun_r <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/phi[2])^2)
}

test_that("UFun is the same in R and TMB", {
  ntest <- 20
  nphi <- 3
  for(ii in 1:ntest) {
    # simulate data
    n <- sample(10:20,1)
    f <- sim_f(n)
    # create TMB model and functions
    tmod <- TMB::MakeADFun(data = list(model_name = "SHOWFit", f = matrix(f)),
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
