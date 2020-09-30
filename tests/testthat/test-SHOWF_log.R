# source("realPSD-testfunctions.R")

# context("SHOWF log parameterization")

# test_that("The SHOWF UFun returned is the same in R and TMB", {
#   ntest <- 20
#   nphi <- sample(2:5,1)
#   for(ii in 1:ntest) {
#     # pick model
#     model <- "SHOWF_log"
#     ufun_r <- get_ufun(model)
#     # simulate data
#     N <- sample(10:20,1)
#     f <- matrix(sim_f(N))
#     phi <- matrix(c(rexp(4), log(0.5+runif(1)) ))
#     # create TMB model and functions
#     tmod <- TMB::MakeADFun(data = list(model = model, method = "UFun", f = f),
#                            parameters = list(phi = phi),
#                            silent = TRUE, DLL = "realPSD_TMBExports")
#     ufun_tmb <- function(phi) c(tmod$simulate(phi)$U)
#     # check they are equal
#     Phi <- replicate(nphi, sim_showf_phi(model = model))
#     U_r <- apply(Phi, 2, ufun_r, f = f)
#     U_tmb <- apply(Phi, 2, ufun_tmb)
#     expect_equal(U_r, U_tmb)
#   }
# })
