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
#' - rename `phi` to `psi`.


#+ gradient_methods

#' Base class for PSD models.
#'
psd_model <- R6::R6Class(
  classname = "psd_mod",

  private = list(

    n_phi_ = NULL, # Number of parameters (excluding scale).
    freq_ = NULL, # Frequency vector.
    Ypsd_ = NULL, # PSD vector.
    fbar_ = NULL, # Binned frequencies.
    Ybar_ = NULL # Binned periodogram ordinates.
    Zbar_ = NULL, # Log binned periodogram.
    Ypsd_scale_ = NULL, # Scaling factor for periodogram.
    Ybar_scale_ = NULL, # Scaling factor for binned periodogram.
    Zbar_scale_ = NULL, # Scaling factor for log binned periodogram.

    # TMB variables
    tmb_DLL_ = NULL, # Path to shared object.
    tmb_model_ = NULL, # Name of model.
    ufun_args_ = NULL, # Additional arguments to UFun constructor.
    ## ufun_ = NULL, # Normalized PSD function.
    ## nlp_ = NULL, # Negative profile loglikelihood.
    ## nll_ = NULL, # Negative loglikelihood.
    ## resid_ = NULL, # Residuals.
    ## zeta_ = NULL, # Conditional estimate of zeta.

    # estimator variables
    est_type_ = NULL, # Type of estimator: lp, nls, mle.
    bin_size_ = NULL, # Bin size.
    bin_type_ = NULL, # Type of binning: mean, median.

    #` Reset bin computations.
    reset_bar = function() {
      sapply(c("fbar_", "Ybar_", "Zbar_"),
             assign, value = NULL, pos = private)
      ## private$fbar_ <- NULL
      ## private$Ybar_ <- NULL
      ## private$Zbar_ <- NULL
    },

    ## #` Reset TMB methods.
    ## reset_tmb = function() {
    ##   sapply(c("ufun_", "nlp_", "nll_", "resid_", "zeta_"),
    ##          assign, value = NULL, pos = private)
    ## },

    #` Prepare TMB `data` list.
    #`
    #` @param method Optional name of the method: "nlp", "nll", "ufun", or "resid".  Gets appended to `est_type`, such that `data$method` becomes e.g., "LP_nlp", "MLE_ufun", etc.
    #` @param ... Additional arguments in `data` list.
    #`
    #` @return A list containing the `data` argument to [TMB::MakeADFun()].
    #`
    #` @note Must be called after `set_est()`.
    get_data = function(method, ...) {
      if(private$est_type_ == "lp") {
        data <- list(
            method = "LP",
            fbar = as.matrix(private$fbar_),
            Zbar = as.matrix(private$Zbar_),
            fs = exp(-bin_factor(private$bin_size_))
        )
      } else if(private$est_type_ == "nls") {
        data <- list(
          method = "NLS",
          fbar = as.matrix(private$fbar_),
          Ybar = as.matrix(private$Ybar_),
          fs = 1
        )
      } else if(private$est_type_ == "mle") {
        data <- list(
            method = "MLE",
            f = as.matrix(private$freq_),
            Y = as.matrix(private$Ypsd_),
            fs = 1
          )
      }
      if(!missing(method)) data$method <- paste0(data$method, "_", method)
      data$model <- private$tmb_model
      data <- c(data, ufun_args_, list(...))
      return(data)
    },

    #` Prepare TMB `parameters` list.
    #`
    #` @param method Name of the method: "nlp", "nll", "ufun", or "resid".      #` @param psi0 Vector with which to initialize `psi`.  If missing, a vector of zeros.
    #` @param zeta0 Scalar with which to initialize `zeta`.  If missing, zero.
    #`
    #` @return A list containing the `parameters` argument to [TMB::MakeADFun()].
    get_parameters = function(method, psi0, zeta0) {
      if(missing(psi0)) psi0 <- rep(0, private$n_phi_)
      if(missing(zeta0)) zeta0 <- 0
      parameters <- list(phi = as.matrix(psi0))
      if(method %in% c("nll", "resid")) {
        # add zeta
        parameters <- c(parameters, list(zeta = zeta0))
        ## if(private$est_type_ == "lp") {
        ##   parameters <- c(parameters, list(zeta = zeta0))
        ## } else if(private$est_type_ == "nls") {
        ##   parameters <- c(parameters, list(tau = exp(zeta0)))
        ## } else if(private$est_type_ == "mle") {
        ##   parameters <- c(parameters, list(tau = exp(zeta0)))
        ## }
      }
      parameters
    }

  ),

  active = list(

    #' @field est_type Estimator type: "lp", "nls", or "mle".
    est_type = function() private$est_type_,

    #' @field bin_size Bin size.
    bin_size = function(value) private$bin_size_,

    #' @field freq Frequency vector.
    freq = function(value) private$freq_,

    #' @field Ypsd Periodogram vector.
    Ypsd = function(value) {
      private$Ypsd_scale_ * private$Ypsd_
    },

    #' @field fbar Binned frequency vector.
    fbar = function(value) {
      private$fbar_
    },

    #' @field Ybar Binned periodogram vector.
    Ybar = function(value) {
      if(is.null(private$Ybar_)) return(NULL)
      private$Ybar_scale_ * private$Ybar_
    },

    #' @field Zbar Log binned periodogram vector.
    Zbar = function(value) {
      if(is.null(private$Zbar_)) return(NULL)
      private$Zbar_scale_ + private$Zbar_
    },

    #' @field n_phi Dimension of `phi`.
    n_phi = function(value) private$n_phi

  ),

  public = list(

    #' @description Specify the estimator.
    #'
    #' @param est_type Type of estimator: "lp", "nls", or "mle".
    #' @param bin_size Bin size.
    #' @param bin_type Type of bin: "mean" or "median".
    set_est = function(est_type = c("lp", "nls", "mle"),
                       bin_size = 100,
                       bin_type = c("mean", "median")) {
      est_type <- match.arg(est_type)
      private$est_type_ <- est_type
      bin_size <- as.integer(bin_size)
      if(length(bin_size) != 1 || bin_size < 1) {
        stop("bin_size must be a positive integer.")
      }
      bin_size <- as.numeric(bin_size)
      private$bin_size_ <- bin_size
      bin_type <- match.arg(bin_type)
      private$bin_type_ <- bin_type
    },

    #' @description Set the periodogram data.
    #'
    #' @param freq Frequency vector.
    #' @param Ypsd Periodogram vector.
    #'
    #' @note Calculates the normalizing constant and binning, so must be called after `set_est()`.
    set_psd = function(freq, Ypsd) {
      N <- length(freq)
      if(length(Ypsd) != N) {
        # todo: better errors
        stop("freq and Ypsd must have same length.")
      }
      private$freq_ <- freq
      private$Ypsd_ <- Ypsd
      # normalize
      private$Ypsd_scale_ <- mean(private$Ypsd_)
      private$Ypsd_ <- private$Ypsd_ / private$Ypsd_scale_
      if(private$est_type %in% c("lp", "nls")) {
        # binning
        private$fbar_ <- binning(private$fseq_,
                                 bin_size = bin_size,
                                 bin_type = bin_type)
        private$Ybar_ <- binning(private$Ypsd_,
                                 bin_size = bin_size,
                                 bin_type = bin_type)
        private$Zbar_ <- log(private$Ybar_)
        # normalize
        private$Ybar_scale <- mean(private$Ybar_)
        private$Zbar_scale <- mean(private$Zbar_)
        private$Ybar_ <- private$Ybar_ / private$Ybar_scale_
        private$Zbar_ <- private$Zbar_ - private$Zbar_scale_
      }
    },

    #' @description Negative profile loglikelihood constructor.
    #'
    #' @param psi0 Vector with which to initialize `psi`.  If missing, a vector of zeros.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `nlp$fn(psi)`, `nlp$gr()`, and `nlp$he()`.
    nlp = function(psi0) {
      method <- "nlp"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, psi0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     silent = TRUE,
                     DLL = private$tmb_DLL)
    },


    #' @description Negative loglikelihood constructor.
    #'
    #' @param psi0 Vector with which to initialize `psi`.  If missing, a vector of zeros.
    #' @param zeta0 Scalar with which to initialize `zeta`.  If missing, zero.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `nll$fn(psi)`, `nll$gr()`, and `nll$he()`.
    nll = function(psi0, zeta0) {
      method <- "nll"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, psi0, zeta0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     silent = TRUE,
                     DLL = private$tmb_DLL)
    },


    #' @description Residual constructor.
    #'
    #' @param psi0 Vector with which to initialize `psi`.  If missing, a vector of zeros.
    #' @param zeta0 Scalar with which to initialize `zeta`.  If missing, zero.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `resid$adreport(psi)$resid` and `resid$gr(psi)`.
    resid = function(psi0, zeta0) {
      method <- "resid"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, psi0, zeta0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     ADreport = TRUE,
                     silent = TRUE,
                     DLL = private$tmb_DLL)
    },

    #' @description Normalized PSD constructor.
    #'
    #' @param psi0 Vector with which to initialize `psi`.  If missing, a vector of zeros.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `ufun$adreport(psi)$U` and `ufun$gr(psi)`.
    ufun = function(psi0) {
      method <- "ufun"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, psi0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     ADreport = TRUE,
                     silent = TRUE,
                     DLL = private$tmb_DLL)
    },

    #' @description Calculate the optimal value of the log normalizing constant.
    #'
    #' @param psi0 Vector with which to initialize `psi`.  If missing, a vector of zeros.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `zeta$fn(psi)`, `zeta$gr()`, and `zeta$he()`.
    zeta = function(psi0) {
      method <- "nlp"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, psi0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     silent = TRUE,
                     DLL = private$tmb_DLL)
    },

    #' @description Calculate the variance estimator for the parameter fitting method.
    #'
    #' @param psi Vector of normalized PSD parameters.
    #' @param zeta Log of normalization parameter.
    #' @param to_theta If `TRUE`, calculates the variance in the inferential basis.  Otherwise, uses the computational basis.
    #'
    #' @return A variance matrix of size `n_phi x n_phi`.
    vcov = function(psi, zeta, to_theta = FALSE) {
      if(private$est_type_ == "lp") {
        nll <- self$nll()
        he <- nll$he(c(psi, zeta))
        out <- 2/private$bin_size_ * chol2inv(chol(he))
      } else if(private$est_type_ == "mle") {
        nll <- self$nll()
        he <- nll$he(c(psi, zeta))
        out <- chol2inv(chol(he))
      } else if(private$est_type_ == "nls") {
        nll <- self$nll()
        resid <- self$resid()
        he <- nll$he(c(psi, zeta))
        ihe <- chol2inv(chol(he))
        jac <- resid$fn(c(psi, zeta)) * resid$grad(c(psi, zeta))
        out <- ihe %*% crossprod(jac) %*% ihe
      }
      if(to_theta) {
        # change of variables
        jac_trans <- numDeriv::jacobian(func = self$itrans,
                                        x = c(psi, zeta))
        out <- jac_trans %*% out %*% t(jac_trans)
      }
      out
    },

    #' @description Model object constructor.
    #'
    #' @param freq Frequency vector.
    #' @param Ypsd Periodogram vector.
    #' @param est_type Type of estimator: "lp", "nls", or "mle".
    #' @param bin_size Bin size.
    #' @param bin_type Type of bin: "mean" or "median".
    #' @param ufun_args Optional list of additional arguments to TMB `UFun()` constructor.
    initialize = function(freq, Ypsd,
                          est_type = c("lp", "nls", "mle"),
                          bin_size = 100,
                          bin_type = c("mean", "median"),
                          ufun_args = NULL) {
      est_type <- match.arg(est_type)
      bin_type <- match.arg(bin_type)
      self$set_est(est_type, bin_size, bin_type)
      self$set_psd(freq, Ypsd)
      private$ufun_args_ <- ufun_args
    }

  )

)
