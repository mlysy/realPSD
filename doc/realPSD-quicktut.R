## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE, 
  message = FALSE,
  fig.align = "center",
  fig.pos = "!htb",
  fig.retina = 2,
  # dpi = 200,
  fig.width = 10, 
  fig.height = 5, 
  out.width = "90%",
  collapse = TRUE,
  comment = "#>"
)
# link to packages
pkg_link <- function(pkg, link) {
  if(link == "github") {
    link <- paste0("https://github.com/mlysy/", pkg)
  } else if(link == "cran") {
    link <- paste0("https://CRAN.R-project.org/package=", pkg)
  }
  paste0("[**", pkg, "**](", link, ")")
}
cran_link <- function(pkg) pkg_link(pkg, "cran")
github_link <- function(pkg) pkg_link(pkg, "github")
# load required packages
require(realPSD)
require(TMB)
require(R6)
require(dplyr)
# get current working directory
wd <- getwd()
# invisibly copy OU_Model.hpp from realPSD to the current working directory
# this trick helps realPSD pass win_builder test in rebuilding vignettes section
file.copy(from = system.file("include", "realPSD", "OU_Model.hpp", 
                              package = "realPSD"),
          to = wd)

## ----ou_psd-------------------------------------------------------------------
#' PSD for the OU model.
#'
#' @param freq Vector of frequencies
#' @param alpha Mean reversion parameter.
#' @param beta Scale parameter.
#' @return PSD vector at the frequencies in `freq`.
ou_psd <- function(freq, alpha, beta) {
  beta^2/((2*pi*freq)^2 + alpha^2)
}

## ----sim----------------------------------------------------------------------
# true parameter values
theta0 <- c(alpha = 3.1, beta = 2.2)

# timescale for simulation
frange <- c(1e-4, 1e2) # frequency range
fs <- frange[2]*2 # sampling frequency
dt <- 1/fs # interobservation time
N <- floor(fs/frange[1]) # number of observations

# time domain simulation
set.seed(314159) # just for reproducibility
Xt <- ou_sim(gamma = theta0[1], mu = 0, sigma = theta0[2], dt = dt, n_obs = N)

## ----psd_emp------------------------------------------------------------------
psd_emp <- periodogram(Xt, fs = fs)
psd_emp <- psd_emp %>% filter(freq < 20)
bin_size <- 100 # for estimation and plotting
psd_bin <- as.data.frame(lapply(psd_emp, binning,
                                bin_size = bin_size, bin_type = "mean"))

plot(psd_bin$freq, psd_bin$Ypsd, type = "l", log = "xy",
     xlab = "Frequency (Hz)", ylab = "PSD")
curve(ou_psd(freq = x, alpha = theta0[1], beta = theta0[2]),
      col = "red", add = TRUE)

## ----ou_cpp_compile_step1-----------------------------------------------------
# create c++ file
ou_cpp <- make_psd_model(model = "OU",
                         header = "OU_Model.hpp",
                         class = "ou::UFun",
                         ctor = "ou::make_Ufun",
                         standalone = TRUE)
# write the content into a .cpp file called OU_FitMethods.cpp under your current working directory
cat(ou_cpp, sep = "\n", file = "OU_FitMethods.cpp")

## ----ou_cpp_compile_step2, message = FALSE------------------------------------
# compile it
model <- "OU_FitMethods"
TMB::compile(paste0(model, ".cpp"),
             PKG_CXXFLAGS = paste0("-I", system.file("include", package = "realPSD"),
                                   " -std=gnu++11"))
dyn.load(TMB::dynlib(model))

## ----ou_cpp_compile_step3-----------------------------------------------------
# create the ou_model class
ou_model <- R6::R6Class(
  classname = "ou_model",
  inherit = realPSD::psd_model,
  private = list(
    n_phi_ = 1,
    tmb_model_ = "OU",
    tmb_DLL_ = "OU_FitMethods"
  ),
  public = list(
    #  Method to convert a value of `phi` into `theta` for the OU model
    #' @param theta Value of `phi = log(alpha)`
    #' @param zeta Value of `zeta = 2*log(beta)`
    #' @return Parameter vector of `theta = c(alpha, beta)`
    to_theta = function(phi, zeta) {
      c(alpha = exp(phi), beta = exp(zeta/2))
    },
    # Method to convert a value of natural parameter `theta` into the log transformed parameter `phi_zeta` for the OU model
    #' @param theta Value of `theta = c(alpha, beta)`
    #' @param obj Model object created by `ou_model`
    #' @return List with element `phi` and `zeta`
    to_phi = function(theta) {
      phi <- log(theta[1])
      zeta <- 2*log(theta[2])
      list(phi = phi, zeta = zeta)
    }
  )
)

## ----tmb_test-----------------------------------------------------------------
# construct the model object
ou_obj <- ou_model$new(freq = psd_emp$freq, Ypsd = psd_emp$Ypsd,
                       bin_size = 100, est_type = "lp")

ou_ufun <- ou_obj$ufun() # TMB method corresponding to UFun.eval()
U_tmb <- ou_ufun$fn(log(theta0[1])) # evaluate at frequencies for lp estimator
# same calculation in R
U_r <- ou_psd(freq = psd_bin$freq, alpha = theta0[1], beta = 1)
range(U_tmb - U_r)

## ----helper_for_debug, echo = FALSE, eval = FALSE-----------------------------
#  #` Helper function to convert a value of `phi_est` into `eta_coef` and `eta_vcov`.
#  #` @param phi Value of `phi`.
#  #` @param obj Model object created by `ou_model`.
#  #` @return List with elements `coef` and `vcov` of parameters on the computational basis, `eta = (phi, zeta)`.
#  to_fit <- function(phi, obj) {
#    fit <- list(phi = phi, zeta = obj$nlp()$simulate(phi)$zeta)
#    fit <- list(coef = c(phi = fit$phi, zeta = fit$zeta),
#                vcov = obj$vcov(phi = fit$phi, zeta = fit$zeta))
#    rownames(fit$vcov) <- colnames(fit$vcov) <- names(fit$coef)
#    fit
#  }
#  # Helper function to convert the fitted parameter `phi` which is on the computational basis (log-transformed) into the natural parameter `theta` for the OU model
#  #' @param phi Value of `phi`
#  #' @param obj Model object created by `ou_model`
#  #' @return A list with elements `coef`, `vcov`, and `se`.
#  to_est <- function(phi, obj) {
#    zeta <- obj$nlp()$simulate(phi)$zeta
#    coef <- obj$to_theta(phi,zeta)
#    vcov <- obj$vcov(phi, zeta, to_theta = TRUE)
#    # change of variables
#    # he <- obj$nll()$he(c(phi, zeta))
#    # out <- chol2inv(chol(he))
#    # jac_trans <- numDeriv::jacobian(func = function(eta) {
#    #         obj$to_theta(phi = eta[1],
#    #                       zeta = eta[2])
#    # }, x = c(phi, zeta))
#    # vcov <- jac_trans %*% out %*% t(jac_trans)
#    colnames(vcov) <- rownames(vcov) <- names(coef)
#    # se <- sqrt(diag(vcov))
#    # list(coef = coef, se = se, vcov = vcov)
#    list(coef = coef, vcov = vcov)
#  }

## ----mle_fit------------------------------------------------------------------
# instantiate TMB object for profile likelihood
ou_obj$set_est(est_type = "mle")
ou_nlp <- ou_obj$nlp()

opt <- optimize(f = ou_nlp$fn, lower = log(.1), upper = log(10))
mle_fit <- ou_obj$to_est(opt$min, to_theta = FALSE) # estimates in the computational/log-transformed basis
print(mle_fit)
mle_est <- ou_obj$to_est(opt$min, to_theta = TRUE) # estimates in the "inferential basis"
print(mle_est)

plot(psd_bin$freq, psd_bin$Ypsd, type = "l", log = "xy",
     xlab = "Frequency (Hz)", ylab = "PSD")
curve(ou_psd(freq = x, alpha = exp(mle_fit$coef[1]), beta = exp(mle_fit$coef[2]/2)),
      col = "red", add = TRUE)

## ----lp_fit-------------------------------------------------------------------
# instantiate TMB object for profile likelihood
ou_obj$set_est(est_type = "lp", bin_size = bin_size)
ou_nlp <- ou_obj$nlp()

opt <- optimize(f = ou_nlp$fn, lower = log(.1), upper = log(10))
lp_fit <- ou_obj$to_est(opt$min) # to_theta = FALSE by default
lp_est <- ou_obj$to_est(opt$min, to_theta = TRUE)

plot(psd_bin$freq, psd_bin$Ypsd, type = "l", log = "xy",
     xlab = "Frequency (Hz)", ylab = "PSD")
curve(ou_psd(freq = x, alpha = exp(lp_fit$coef[1]), beta = exp(lp_fit$coef[2]/2)),
      col = "red", add = TRUE)

## ----nls_fit------------------------------------------------------------------
# instantiate TMB object for profile likelihood
ou_obj$set_est(est_type = "nls", bin_size = bin_size)
ou_nlp <- ou_obj$nlp()

opt <- optimize(f = ou_nlp$fn, lower = log(.1), upper = log(10))
nls_fit <- ou_obj$to_est(opt$min)
nls_est <- ou_obj$to_est(opt$min, to_theta = TRUE)

plot(psd_bin$freq, psd_bin$Ypsd, type = "l", log = "xy",
     xlab = "Frequency (Hz)", ylab = "PSD")
curve(ou_psd(freq = x, alpha = exp(lp_fit$coef[1]), beta = exp(lp_fit$coef[2]/2)),
      col = "red", add = TRUE)

## ----disp---------------------------------------------------------------------
disp <- sapply(list(mle = mle_fit, lp = lp_fit, nls = nls_fit),
               function(fit) c(fit$coef[1], sqrt(fit$vcov[1,1]),
                               fit$coef[2], sqrt(fit$vcov[2,2])))
disp <- cbind(true = c(log(theta0[1]), NA, 2 * log(theta0[2]), NA),
              disp)
disp <- t(disp)
colnames(disp) <- c("phi_est", "phi_se", "zeta_est", "zeta_se")
signif(disp, 3)

## ----disp2--------------------------------------------------------------------
disp <- sapply(list(mle = mle_est, lp = lp_est, nls = nls_est),
               function(fit) c(fit$coef[1], sqrt(fit$vcov[1,1]),
                               fit$coef[2], sqrt(fit$vcov[2,2])))
disp <- cbind(true = c(theta0[1], NA, theta0[2], NA), disp)
disp <- t(disp)
colnames(disp) <- c("alpha_est", "alpha_se", "beta_est", "beta_se")
signif(disp, 3)

## ----ou_ufun, echo = FALSE, results = "asis"----------------------------------
cat("```cpp", readLines("OU_Model.hpp"), "```", sep = "\n")

