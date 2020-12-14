#' ---
#' title: "Application: Calibration of an Atomic Force Microscope"
#' ---

# packages

require(tidyverse)
require(realPSD)
require(optimCheck)
require(numDeriv)
require(R6)

#' ## data and psd

#+ data
psd_data <- readRDS("data/psd_data.rds")
SF <- 5e6 # sampling frequency 5MHz
T_s <- 5 # total time (second)
Temp <- 298
scale_factor <- 1
psd_data <- periodogram(psd_data$ytime, SF_s = SF, T_s = T_s)
psd_data <- psd_data %>% as_tibble() %>%
  mutate(yFreq = yFreq/(2*pi*xFreq / 0.05)^2) %>% # normalization
  mutate(yFreq = yFreq/2) # unfolding
xPSD <- psd_data$xFreq
yPSD <- psd_data$yFreq * scale_factor

#' ## fit NLS estimator in pure R

#+ nls_r

# Calculate SHOF PSD.
shof_psd <- function(fseq, log_theta, Temp) {
  theta <- as.list(exp(log_theta))
  ## names(theta)[names(theta) == "Sf"] <- "Af"
  do.call(showf_psd,
          c(list(fseq = fseq, Temp = Temp, Sw = 0), theta))
  ## showf_psd(fseq, k, f0, Q, 0, Sf, alpha, Temp)
}

# initial values
theta0 <- c(k = .178/scale_factor, f0 = 33563.03, Q = 52.8,
            Af = 7.218499e-18 * scale_factor,
            alpha = 1.99)
log_theta0 <- log(theta0)

# frequency range
frange <- c(1, 60)*1000
findex <- c(which.min(abs(frange[1] - xPSD)) : which.min(abs(frange[2] - xPSD)))

# denoising
yPSD_denoise <- psd_denoise(Yf = yPSD[findex],
                            psd_fs = shof_psd(fseq = xPSD[findex], Temp = Temp,
                                              log_theta = log_theta0),
                            method = "fisherG")


#' ### sanity check

bin_size <- 100
fbar <- binning(xPSD[findex], bin_size = bin_size, bin_type = "mean")
## Ybar <- binning(yPSD[findex], bin_size = bin_size, bin_type = "mean")
Ybar <- binning(yPSD_denoise, bin_size = bin_size, bin_type = "mean")

plot(fbar, Ybar, type = "l", log = "xy")
lines(fbar,
      shof_psd(fseq = fbar, Temp = Temp, log_theta = log_theta0),
      col = "red")

#' ### now fit NLS model

resid_fun <- function(log_theta) {
  constY <- mean(Ybar)^2
  (Ybar - shof_psd(fbar, log_theta, Temp))^2/constY
}
obj_fun <- function(log_theta) sum(resid_fun(log_theta))

fit <- optim(par = log_theta0, fn = obj_fun,
             control = list(maxit = 10000
                            ## reltol = 1e-20,
                            ## parscale = theta0
                            ))

#' ### sanity check

# all fine except alpha
optim_proj(xsol = fit$par,
           fun = obj_fun, maximize = FALSE)

#' ### sandwich estimator

# variance of log_theta_hat
he <- hessian(func = obj_fun, x = fit$par)
ihe <- chol2inv(chol(he))
jac <- jacobian(func = resid_fun, x = fit$par)
vcov <- ihe %*% crossprod(jac) %*% ihe

# variance of theta_hat: change of variables
J <- jacobian(func = exp, x = fit$par)
sqrt(diag(J %*% vcov %*% J))

#' ## Estimators and Standard Errors
#'
#' In order to use TMB and R effectively, we propose:
#'
#' - Likelihood gradients with TMB in a favored computational basis (e.g., unconstrained in $\mathbb{R}^d$).
#' - Numerical gradients between computational and inferential bases in R using **numDeriv**.
#' - As methods of a class can implement TMB gradients for these.
#'
#'
#' ### Notation
#'
#' - Inferential basis: $\tth = (\pph, \sigma)$.
#' - Computational basis: $\oom = (\pps, \zeta = log \sigma^2)$.
#'
#' ### Proposed changes to TMB functions
#'
#' These changes can be done in `FitMethods.hpp` to avoid making changes all over the place.
#'
#' - `nlp`, `nll`, `resid`, `zeta` functions should have a scale factor which allows them to compute the actual functions (to avoid confusion).
#' - `resid` and `nll` methods should always take `zeta` as the normalizing constant.
#' - conditional optimum should always be `zeta(psi)`.
#' - add residual function for MLE.  Let's say its `log(Ypsd/(tau * Ufun(psi)))`.
#' - add `ufun` method to estimators.
#' - remove `fs` argument from TMB interface.
#' - rename `phi` to `psi`. actually, don't do this.  `theta = (phi, zeta = log sigma^2)` already suggests that `(phi, zeta)` is the computational basis.
#'
#' ### Testing
#'
#' Adapt the existing tests.  That is, the TMB parts should be created by the R6 class.


#+ gradient_methods

