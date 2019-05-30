sim_f <- function(n) runif(n, 0, 2*n)
sim_Zbar <- function(n) rnorm(n)
sim_phi <- function() c(f0 = runif(1, 100, 1000), gamma = rexp(1), Rw = rexp(1))
sim_zeta <- function() rexp(1)
zeta_r <- function(fbar, Zbar, phi) {
  mean(Zbar - log(ufun_r(fbar, phi)))
}
ufun_r <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/phi[2])^2)
}
nllik_r <- function(phi, zeta, Zbar, fbar) {
  logUbar <- log(ufun_r(fbar, phi))
  sum((Zbar - zeta - logUbar)^2)
}
tmb_recompile <- function() {
  RcppTMBTest::export_models()
  pkgbuild::compile_dll()
  devtools::install(quick = TRUE)
  q()
}
