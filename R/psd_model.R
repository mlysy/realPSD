#' Base class for PSD models.
#'
#' @export
psd_model <- R6::R6Class(
  classname = "psd_model",

  private = list(

    n_phi_ = NULL, # Number of parameters (excluding scale).
    freq_ = NULL, # Frequency vector.
    Ypsd_ = NULL, # PSD vector.
    fbar_ = NULL, # Binned frequencies.
    Ybar_ = NULL, # Binned periodogram ordinates.
    Zbar_ = NULL, # Log binned periodogram.
    Ypsd_scale_ = NULL, # Scaling factor for periodogram.
    Ybar_scale_ = NULL, # Scaling factor for binned periodogram.
    Zbar_scale_ = NULL, # Scaling factor for log binned periodogram.

    # TMB variables
    tmb_DLL_ = NULL, # Path to shared object.
    tmb_model_ = NULL, # Name of model.
    ctor_args_ = NULL, # Additional arguments to UFun constructor.
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

    #` Calculate bin factors for LP and NLS estimators.
    #`
    #` The bin factor is the expected value of the bin variable (regular/log mean/median), and in the case of LP, also the `nll` scaling factor.
    get_scale = function() {
      est_type <- private$est_type_
      bin_size <- private$bin_size_
      bin_type <- private$bin_type_
      if(est_type == "lp") {
        scale <- private$Zbar_scale_
        if(bin_type == "mean") {
          scale <- c(scale - (digamma(bin_size) - log(bin_size)),
                  (bin_size/2) / (bin_size * trigamma(bin_size)))
        } else if(bin_type == "median") {
          scale <- c(scale - log(log(2)), 1)
        }
      } else if(est_type == "nls") {
        scale <- log(private$Ybar_scale_)
        if(bin_type == "median") {
          scale <- scale - log(log(2))
        }
      } else if(est_type == "mle") {
        scale <- log(private$Ypsd_scale_)
      }
      scale
    },

    #` Prepare TMB `data` list.
    #`
    #` @param method Optional name of the method: "nlp", "nll", "ufun", or "resid".  Gets appended to `est_type`, such that `data$method` becomes e.g., "LP_nlp", "MLE_ufun", etc.
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
            scale = private$get_scale()
        )
      } else if(private$est_type_ == "nls") {
        data <- list(
          method = "NLS",
          fbar = as.matrix(private$fbar_),
          Ybar = as.matrix(private$Ybar_),
          scale = private$get_scale()
        )
      } else if(private$est_type_ == "mle") {
        data <- list(
            method = "MLE",
            f = as.matrix(private$freq_),
            Y = as.matrix(private$Ypsd_),
            scale = private$get_scale()
          )
      }
      if(!missing(method)) data$method <- paste0(data$method, "_", method)
      data$model <- private$tmb_model_
      data <- c(data, private$ctor_args_)
      return(data)
    },

    #` Prepare TMB `parameters` list.
    #`
    #` @param method Name of the method: "nlp", "nll", "ufun", or "resid".      #` @param phi0 Vector with which to initialize `phi`.  If missing, a vector of zeros.
    #` @param zeta0 Scalar with which to initialize `zeta`.  If missing, zero.
    #`
    #` @return A list containing the `parameters` argument to [TMB::MakeADFun()].
    get_parameters = function(method, phi0, zeta0) {
      if(missing(phi0)) phi0 <- rep(0, private$n_phi_)
      if(missing(zeta0)) zeta0 <- 0
      parameters <- list(phi = as.matrix(phi0))
      if(method %in% c("nll", "res")) {
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
      if(private$est_type_ %in% c("lp", "nls")) {
        # binning
        private$fbar_ <- binning(private$freq_,
                                 bin_size = private$bin_size_,
                                 bin_type = private$bin_type_)
        private$Ybar_ <- binning(private$Ypsd_,
                                 bin_size = private$bin_size_,
                                 bin_type = private$bin_type_)
        private$Zbar_ <- log(private$Ybar_)
      }
      # normalize
      private$Ypsd_scale_ <- mean(private$Ypsd_)
      private$Ypsd_ <- private$Ypsd_ / private$Ypsd_scale_
      if(private$est_type_ %in% c("lp", "nls")) {
        private$Ybar_scale_ <- mean(private$Ybar_)
        private$Zbar_scale_ <- mean(private$Zbar_)
        private$Ybar_ <- private$Ybar_ / private$Ybar_scale_
        private$Zbar_ <- private$Zbar_ - private$Zbar_scale_
      }
    },

    #' @description Negative profile loglikelihood constructor.
    #'
    #' @param phi0 Vector with which to initialize `phi`.  If missing, a vector of zeros.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `nlp$fn(phi)`, `nlp$gr(phi)`, `nlp$he(phi)`, and `nlp$simulate(phi)$zeta`, which returns the conditional optimum of `zeta` given `phi`.
    nlp = function(phi0) {
      method <- "nlp"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, phi0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     silent = TRUE,
                     DLL = private$tmb_DLL_)
    },


    #' @description Negative loglikelihood constructor.
    #'
    #' @param phi0 Vector with which to initialize `phi`.  If missing, a vector of zeros.
    #' @param zeta0 Scalar with which to initialize `zeta`.  If missing, zero.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `nll$fn(phi)`, `nll$gr()`, and `nll$he()`.
    nll = function(phi0, zeta0) {
      method <- "nll"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, phi0, zeta0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     silent = TRUE,
                     DLL = private$tmb_DLL_)
    },


    #' @description Residual constructor.
    #'
    #' @param phi0 Vector with which to initialize `phi`.  If missing, a vector of zeros.
    #' @param zeta0 Scalar with which to initialize `zeta`.  If missing, zero.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `resid$fn(phi)` and `resid$gr(phi)`.
    resid = function(phi0, zeta0) {
      method <- "res"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, phi0, zeta0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     ADreport = TRUE,
                     silent = TRUE,
                     DLL = private$tmb_DLL_)
    },

    #' @description Normalized PSD constructor.
    #'
    #' @param phi0 Vector with which to initialize `phi`.  If missing, a vector of zeros.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `ufun$fn(phi)` and `ufun$gr(phi)`.
    ufun = function(phi0) {
      method <- "ufun"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, phi0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     ADreport = TRUE,
                     silent = TRUE,
                     DLL = "OU_FitMethods")
    },

    #' @description Calculate the optimal value of the log normalizing constant.
    #'
    #' @param phi0 Vector with which to initialize `phi`.  If missing, a vector of zeros.
    #'
    #' @return A [TMB::MakeADFun()] object.  Relevant methods are: `zeta$fn(phi)`, `zeta$gr()`, and `zeta$he()`.
    zeta = function(phi0) {
      method <- "zeta"
      data <- private$get_data(method)
      parameters <- private$get_parameters(method, phi0)
      TMB::MakeADFun(data = data,
                     parameters = parameters,
                     silent = TRUE,
                     DLL = private$tmb_DLL_)
    },

    #' @description Calculate the variance estimator for the parameter fitting method.
    #'
    #' @param phi Vector of normalized PSD parameters.
    #' @param zeta Log of normalization parameter.
    #' @param to_theta If `TRUE`, calculates the variance in the inferential basis.  Otherwise, uses the computational basis.
    #'
    #' @return A variance matrix of size `n_phi x n_phi`.
    vcov = function(phi, zeta, to_theta = FALSE) {
      if(private$est_type_ == "lp") {
        nll <- self$nll()
        he <- nll$he(c(phi, zeta))
        out <- chol2inv(chol(he))
      } else if(private$est_type_ == "mle") {
        nll <- self$nll()
        he <- nll$he(c(phi, zeta))
        out <- chol2inv(chol(he))
      } else if(private$est_type_ == "nls") {
        nll <- self$nll()
        resid <- self$resid()
        he <- nll$he(c(phi, zeta))
        ihe <- chol2inv(chol(he))
        jac <- 2 * resid$fn(c(phi, zeta)) * resid$gr(c(phi, zeta))
        out <- ihe %*% crossprod(jac) %*% ihe
      }
      if(to_theta) {
        # change of variables
        jac_trans <- numDeriv::jacobian(func = self$itrans,
                                        x = c(phi, zeta))
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
    #' @param ctor_args Optional list of additional arguments to TMB `UFun()` constructor.
    initialize = function(freq, Ypsd,
                          est_type = c("lp", "nls", "mle"),
                          bin_size = 100,
                          bin_type = c("mean", "median"),
                          ctor_args = NULL) {
      est_type <- match.arg(est_type)
      bin_type <- match.arg(bin_type)
      self$set_est(est_type, bin_size, bin_type)
      self$set_psd(freq, Ypsd)
      private$ctor_args_ <- ctor_args
    }

  )

)
