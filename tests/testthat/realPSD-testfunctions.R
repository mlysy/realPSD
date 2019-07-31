# simulation functions
sim_f <- function(n) runif(n, 0, 2*n)
sim_Zbar <- function(n) rnorm(n)
sim_Y <- function(n) rchisq(n, df = 2)
sim_phi <- function() c(f0 = runif(1, 100, 1000), gamma = rexp(1), Rw = rexp(1))
sim_zeta <- function() rexp(1)
sim_tau <- function() rexp(1)
sim_fs <- function() sample(10:1000, size = 1)

# lp functions
lp_zeta_r <- function(fbar, Zbar, phi, ufun, fs) {
  logUbar <- log(fs * ufun(fbar, phi))
  mean(Zbar - logUbar)
}
lp_nll_r <- function(phi, zeta, Zbar, fbar, ufun, fs) {
  logUbar <- log(fs * ufun(fbar, phi))
  sum((Zbar - zeta - logUbar)^2)
}
lp_nlp_r <- function(phi, Zbar, fbar, ufun, fs) {
  zeta <- lp_zeta_r(fbar, Zbar, phi, ufun, fs)
  lp_nll_r(phi, zeta, Zbar, fbar, ufun, fs)
}

# mle functions
mle_nll_r <- function(phi, tau, Y, f, ufun) {
  U <- ufun(f, phi)
  S <- tau * U
  sum(Y/S + log(S))
}
mle_tau_r <- function(f, Y, phi, ufun) {
  U <- ufun(f, phi)
  mean(Y/U)
}
mle_nlp_r <- function(phi, Y, f, ufun) {
  tau <- mle_tau_r(f, Y, phi, ufun)
  mle_nll_r(phi, tau, Y, f, ufun)
  ## length(Y) * (1 + log(tau)) + sum(log(U))
}

# nls functions
nls_nll_r <- function(phi, tau, Ybar, fbar, ufun, fs) {
  Ubar <- fs * ufun(fbar, phi)
  sum((Ybar - tau * Ubar)^2)
}
nls_tau_r <- function(fbar, Ybar, phi, ufun, fs) {
  Ubar <- fs * ufun(fbar, phi)
  sum(Ybar * Ubar) / sum(Ubar * Ubar)
}
nls_nlp_r <- function(phi, Ybar, fbar, ufun, fs) {
  tau <- nls_tau_r(fbar, Ybar, phi, ufun, fs)
  nls_nll_r(phi, tau, Ybar, fbar, ufun, fs)
}

# show model normalized PSD
show_ufun <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/phi[2])^2)
}

# recompile TMB models, install package, and quit
tmb_recompile <- function() {
  RcppTMBTest::export_models()
  pkgbuild::compile_dll()
  devtools::install(quick = TRUE)
  q()
}
