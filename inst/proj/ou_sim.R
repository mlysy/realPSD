#' Generate observations from an Ornstein-Uhlenbeck process at regular time intervals.
#'
#' @param gamma Scalar mean reversion parameter (see **Details**).
#' @param mu Scalar mean parameter (see **Details**).
#' @param sigma Scalar diffusion parameter (see **Details**).
#' @param dt Interobservation time.
#' @param n_obs Number of observations to generate.
#' @param X0 Initial value of the process at time `t = 0`.  If missing sampled from the OU stationary distribution.
#' @return A vector of `n_obs` observations of the process at times `t = dt, 2*dt, ..., n_obs * dt`.
#' @details The Ornstein-Uhlenbeck (OU) process satisfies a stochastic differential equation of the form
#' ```
#' dX_t = -gamma * (X_t - mu) dt + sigma dB_t.
#' ```
#' It is a stationary Gaussian Markov process with transition density
#' ```
#' X_s+t | X_s ~ N( rho_t * (X_s - mu) + mu, tau^2 * (1-rho_t^2) ),
#' ```
#' where `rho_t = exp(-gamma * t)` and `tau^2 = sigma^2/(2*gamma)`.  Its stationary distribution is `X_t ~ N(mu, tau^2)`.
#' @export
ou_sim <- function(gamma, mu, sigma, dt, n_obs, X0) {
  tau <- sigma/sqrt(2*gamma) # stationary standard deviation
  if(missing(X0)) X0 <- rnorm(1, mean = mu, sd = tau)
  # generate efficiently using a one-step linear filter
  lrho <- -gamma * dt
  nu <- tau * sqrt((1-exp(2 * lrho))) # conditional sd
  rho <- exp(lrho) # filter coefficients
  z <- rnorm(n_obs, sd = nu) # pre-generate normal draws
  X <- stats::filter(x = z, filter = rho,
                     method = "recursive", init = X0 - mu)
  as.numeric(X + mu)
}
