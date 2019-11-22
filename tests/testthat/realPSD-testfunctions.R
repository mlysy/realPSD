# simulation functions
sim_f <- function(n) runif(n, 0, 2*n)
sim_Zbar <- function(n) rnorm(n)
sim_Y <- function(n) rchisq(n, df = 2)
sim_phi <- function(model = c("SHOW_nat", "SHOW_log")) {
  model <- match.arg(model)
  phi <- c(f0 = runif(1, 100, 1000),
           Q = runif(1, 1, 500),
           Rw = rexp(1))
  if(model == "SHOW_log") phi <- log(phi)
  phi
}
sim_zeta <- function() rexp(1)
sim_tau <- function() rexp(1)
sim_fs <- function() sample(10:1000, size = 1)
sim_model <- function() sample(c("SHOW_nat", "SHOW_log"), size = 1)

# lp functions
lp_zeta_r <- function(fbar, Zbar, phi, ufun, fs) {
  logUbar <- log(fs * ufun(fbar, phi))
  mean(Zbar - logUbar)
}
lp_res_r <- function(fbar, Zbar, phi, ufun, fs) {
  logUbar <- log(fs * ufun(fbar, phi))
  zeta <- lp_zeta_r(fbar, Zbar, phi, ufun, fs)
  Zbar - zeta - logUbar
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
mle_nll_r <- function(phi, tau, Y, f, ufun, fs) {
  U <- fs * ufun(f, phi)
  S <- tau * U
  sum(Y/S + log(S))
}
mle_tau_r <- function(f, Y, phi, ufun, fs) {
  U <- fs * ufun(f, phi)
  mean(Y/U)
}
mle_nlp_r <- function(phi, Y, f, ufun, fs) {
  tau <- mle_tau_r(f, Y, phi, ufun, fs)
  mle_nll_r(phi, tau, Y, f, ufun, fs)
  ## length(Y) * (1 + log(tau)) + sum(log(U))
}

# nls functions
nls_nll_r <- function(phi, tau, Ybar, fbar, ufun, fs) {
  Ubar <- fs * ufun(fbar, phi)
  sum((Ybar - tau * Ubar)^2)
}
nls_res_r <- function(phi, Ybar, fbar, ufun, fs) {
  Ubar <- fs * ufun(fbar, phi)
  tau <- nls_tau_r(fbar, Ybar, phi, ufun, fs)
  Ybar - tau * Ubar
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
# phi = c(f0, Q, Rw)
show_ufun <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/(phi[1]*phi[2]))^2)
}

# pick model
get_ufun <- function(model = c("SHOW_nat", "SHOW_log")) {
  model <- match.arg(model)
  if(model == "SHOW_nat") {
    return(show_ufun)
  } else if(model == "SHOW_log") {
    return(function(f, phi) show_ufun(f, exp(phi)))
  }
}

# recompile TMB models, install package, and quit
tmb_recompile <- function() {
  RcppTMBTest::export_models()
  pkgbuild::compile_dll()
  devtools::install(quick = TRUE)
  q()
}
