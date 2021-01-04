#' SHOF model class.
#'
#' @export
shof_model <- R6::R6Class(
  classname = "shof_model",
  inherit = psd_model,

  private = list(
    Temp_ = NULL, # temperature
    n_phi_ = 4,
    tmb_DLL_ = "realPSD_TMBExports",
    tmb_model_ = "SHOF_log"
  ),

  active = list(
    #' @field Temp Temperature of the system.
    Temp = function() private$Temp_
  ),

  public = list(
    #' @description Transform computational to original basis.
    #'
    #' @param phi Vector of normalized PSD parameters.
    #' @param zeta Log-scale parameter.
    #'
    #' @return The parameter vector in the original basis.
    #' @details The original basis parameters are
    #'
    #' ```
    #' theta = (k, f0, Q, Sf, alpha).
    #' ```
    #' The transformed parameters are
    #' ```
    #' phi = (log(f0), log(Q), log(Rf), alpha),   zeta = log(tau),
    #' ```
    #' where
    #' ```
    #' tau = Kb*Temp / (k*pi*f0*Q),    Rf = Sf/tau.
    #' ```
    to_theta = function(phi, zeta) {
      Kb <- 1.381e-23
      f0 <- exp(phi[1])
      Q <- exp(phi[2])
      tau <- exp(zeta)
      Sf <- exp(phi[3]) * tau
      alpha <- phi[4]
      k <- Kb*private$Temp_ / (tau * pi * f0 * Q)
      setNames(c(k, f0, Q, Sf, alpha),
               nm = c("k", "f0", "Q", "Sf", "alpha"))
    },

    #' @description Transformation from the original to the computational basis.
    #'
    #' @param theta Vector of parameters in the original basis.
    #' @return A list with elements `phi` and `zeta`.
    to_phi = function(theta) {
      Kb <- 1.381e-23
      log_f0 <- log(theta[2])
      log_Q <- log(theta[3])
      tau <- Kb * private$Temp_ / (theta[1] * pi * theta[2] * theta[3])
      zeta <- log(tau)
      log_Rf <- log(theta[4]) - zeta
      alpha <- theta[5]
      phi <- setNames(c(log_f0, log_Q, log_Rf, alpha),
                      nm = c("log_f0", "log_Q", "log_Rf", "alpha"))
      list(phi = phi, zeta = as.numeric(zeta))
    },

    #' @description SHOF model constructor.
    #'
    #' @param freq Frequency vector.
    #' @param Ypsd Periodogram vector.
    #' @param Temp Temperature (Kelvin).
    #' @param est_type Type of estimator: "lp", "nls", or "mle".
    #' @param bin_size Bin size.
    #' @param bin_type Type of bin: "mean" or "median".
    initialize = function(freq, Ypsd, Temp,
                          est_type = c("lp", "nls", "mle"),
                          bin_size = 100,
                          bin_type = c("mean", "median")) {
      private$Temp_ <- Temp
      super$initialize(freq, Ypsd, est_type, bin_size, bin_type)
    }
  )
)
