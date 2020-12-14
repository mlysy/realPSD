# simulation functions
sim_setup <- function(nphi, est_type = c("lp", "nls", "mle"),
                    bin_type = c("mean", "median")) {
  est_type <- match.arg(est_type)
  bin_type <- match.arg(bin_type)
  # pick model
  model <- sim_model()
  ufun_r <- get_ufun(model)
  # simulate data
  N <- sample(10:20,1)
  B <- sample(2:5, 1)
  freq <- sort(sim_f(N))
  Ypsd <- sim_Y(N)
  fbar <- realPSD::binning(freq, bin_size = B, bin_type = bin_type)
  Ybar <- realPSD::binning(Ypsd, bin_size = B, bin_type = bin_type)
  Zbar <- log(Ybar)
  # simulate parameters
  Phi <- replicate(nphi, sim_phi())
  zeta <- replicate(nphi, sim_zeta())
  # create TMB model class
  test_model <-
    R6::R6Class(
          classname = "test_model",
          inherit = realPSD::psd_model,
          private = list(
            Temp_ = NULL, # temperature
            n_phi_ = 3,
            tmb_DLL_ = "realPSD_TMBExports",
            tmb_model_ = model
          ),
          active = list(
            #' @field Temp Temperature of the system.
            Temp = function() private$Temp_
          )
        )
  # calculate bin constants
  if(est_type == "mle") {
    bin_loc <- 1
  } else if(est_type == "lp") {
    bin_loc <- switch(bin_type,
                       mean = digamma(B) - log(B),
                       median = log(log(2)))
  } else if(est_type == "nls") {
    bin_loc <- switch(bin_type,
                       mean = 1,
                       median = log(2))
  }
  if(est_type == "lp" && bin_type == "mean") {
    bin_scale <- B/2 / (B * trigamma(B))
  } else bin_scale <- 1
  list(model = model, ufun_r = ufun_r, N = N, B = B,
       freq = freq, Ypsd = Ypsd, fbar = fbar, Ybar = Ybar, Zbar = Zbar,
       Phi = Phi, zeta = zeta,
       bin_loc = bin_loc, bin_scale = bin_scale,
       test_model = test_model)
}
sim_f <- function(n) runif(n, 0, 2*n)
sim_Zbar <- function(n) rnorm(n)
sim_Y <- function(n) rchisq(n, df = 2)
sim_phi <- function(model = c("SHOW_nat", "SHOW_log", "SHOW_comp")) {
  model <- match.arg(model)
  phi <- c(f0 = runif(1, 100, 1000),
           Q = runif(1, 1, 500),
           Rw = rexp(1))
  if(model == "SHOW_log") phi <- log(phi)
  phi
}
# TODO: add SHOWF_nat model later on
sim_showf_phi <- function(model = c("SHOWF_log", "SHOWF_nat")) {
  model <- match.arg(model)
  phi <- c(
    f0 = runif(1, 100, 1000),
    Q = runif(1, 1, 500),
    Rw = rexp(1),
    Rf = rexp(1),
    alpha = 0.5+runif(1) # usually 0.5 <= alpha <= 1.5
  )
  if(model == "SHOWF_log") phi <- log(phi)
  phi
}
sim_zeta <- function() log(rexp(1))
sim_tau <- function() rexp(1)
sim_fs <- function() sample(10:1000, size = 1)
sim_model <- function() sample(c("SHOW_nat", "SHOW_log", "SHOW_comp"), size = 1)
## sim_model <- function() {"SHOW_comp"} # force the unit test to test SHOW_comp only

# lp functions
lp_zeta_r <- function(fbar, Zbar, phi, ufun, fs, B, bin_type) {
  bin_loc <- switch(bin_type,
                    mean = digamma(B) - log(B),
                    median = log(log(2)))
  logUbar <- log(fs * ufun(fbar, phi))
  mean(Zbar - logUbar) - bin_loc
}
lp_res_r <- function(fbar, Zbar, phi, zeta, ufun, fs, B, bin_type) {
  bin_loc <- switch(bin_type,
                    mean = digamma(B) - log(B),
                    median = log(log(2)))
  logUbar <- log(fs * ufun(fbar, phi))
  ## zeta <- lp_zeta_r(fbar, Zbar, phi, ufun, fs)
  (Zbar - bin_loc) - zeta - logUbar
}
lp_nll_r <- function(phi, zeta, Zbar, fbar, ufun, fs, B, bin_type) {
  bin_loc <- switch(bin_type,
                    mean = digamma(B) - log(B),
                    median = log(log(2)))
  bin_scale <- switch(bin_type,
                      mean = B/2 / (B * trigamma(B)),
                      median = 1)
  logUbar <- log(fs * ufun(fbar, phi))
  bin_scale * sum(((Zbar - bin_loc) - zeta - logUbar)^2)
}
lp_nlp_r <- function(phi, Zbar, fbar, ufun, fs, B, bin_type) {
  zeta <- lp_zeta_r(fbar, Zbar, phi, ufun, fs, B, bin_type)
  lp_nll_r(phi, zeta, Zbar, fbar, ufun, fs, B, bin_type)
}
# lp function, gradient and hessian
lp_zeta_gr_r <- function(phi, fbar, Zbar, ufun, fs, B) {
  numDeriv::grad(func = lp_zeta_r, x = phi,
                  fbar = fbar, Zbar = Zbar,
                  ufun = ufun, fs = fs, B = B)
}
lp_zeta_he_r <- function(phi, fbar, Zbar, ufun, fs, B) {
  numDeriv::hessian(func = lp_zeta_r, x = phi,
                    fbar = fbar, Zbar = Zbar,
                    ufun = ufun, fs = fs, B = B)
}

lp_nll_gr_r <- function(phi, zeta, fbar, Zbar, ufun, fs, B) {
  numDeriv::grad(func = lp_nll_r, x = phi,
                  zeta = zeta, fbar = fbar, Zbar = Zbar,
                  ufun = ufun, fs = fs, B = B)
}
lp_nll_he_r <- function(phi, zeta, fbar, Zbar, ufun, fs, B) {
  numDeriv::hessian(func = lp_nll_r, x = phi,
                    zeta = zeta, fbar = fbar, Zbar = Zbar,
                    ufun = ufun, fs = fs, B = B)
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
nls_nll_r <- function(phi, tau, Ybar, fbar, ufun, fs, bin_type) {
  bin_loc <- switch(bin_type,
                    mean = 1,
                    median = log(2))
  Ubar <- fs * ufun(fbar, phi)
  sum((Ybar - bin_loc * tau * Ubar)^2)
}
nls_res_r <- function(phi, tau, Ybar, fbar, ufun, fs, bin_type) {
  bin_loc <- switch(bin_type,
                    mean = 1,
                    median = log(2))
  Ubar <- fs * ufun(fbar, phi)
  ## tau <- nls_tau_r(fbar, Ybar, phi, ufun, fs)
  Ybar - bin_loc * tau * Ubar
}
nls_tau_r <- function(fbar, Ybar, phi, ufun, fs, bin_type) {
  bin_loc <- switch(bin_type,
                    mean = 1,
                    median = log(2))
  Ubar <- bin_loc * fs * ufun(fbar, phi)
  sum(Ybar * Ubar) / sum(Ubar * Ubar)
}
nls_nlp_r <- function(phi, Ybar, fbar, ufun, fs, bin_type) {
  tau <- nls_tau_r(fbar, Ybar, phi, ufun, fs, bin_type)
  nls_nll_r(phi, tau, Ybar, fbar, ufun, fs, bin_type)
}

# show model normalized PSD
# phi = c(f0, Q, Rw)
show_ufun <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/(phi[1]*phi[2]))^2)
}

# show model normalized PSD with another parametrization
# phi = c(f0, gamma = f0*Q, Rw)
show_ufun_comp <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/phi[2])^2)
}

# showf model with natural parametrization
# phi = c(f0, Q, Rw, Rf, alpha)
showf_ufun <- function(f, phi) {
  phi[3] + 1/(((f/phi[1])^2 - 1)^2 + (f/(phi[1]*phi[2]))^2) + phi[4]/f^phi[5]
}

# pick model
get_ufun <- function(model = c("SHOW_nat", "SHOW_log", "SHOW_comp", "SHOWF_log", "SHOWF_nat")) {
  model <- match.arg(model)
  if(model == "SHOW_nat") {
    return(show_ufun)
  } else if(model == "SHOW_log") {
    return(function(f, phi) show_ufun(f, exp(phi)))
  } else if(model == "SHOW_comp") {
    return(show_ufun_comp)
  } else if(model == "SHOWF_log") {
    return(function(f, phi) showf_ufun(f, exp(phi)))
  } else if(model == "SHOWF_nat") {
    return(showf_ufun)
  }
}

# recompile TMB models, install package, and quit
tmb_recompile <- function() {
  TMBtools::export_models()
  pkgbuild::compile_dll()
  devtools::document()
  devtools::install(quick = TRUE)
  q()
}

tmb_restart <- function() {
  require(realPSD)
  require(testthat)
}
