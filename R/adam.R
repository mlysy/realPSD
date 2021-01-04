#' ADAM optimizer.
#'
#' @param par Initial parameter value (numeric vector).
#' @param fn Objective function to minimize.
#' @param gr Gradient function of `fn`.
#' @param nsteps Number of optimization steps to perform.
#' @param alpha,beta,eps Algorithm tuning parameters.  `alpha` and `eps` are scalars, `beta` is a numeric vector of length 2.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`par`}{The parameter value at the minimum of the `nsteps` optimization steps.}
#'   \item{`objective`}{The minimum value of the objective function.}
#'   \item{`par_full`}{A matrix where each row is the parameter value at the given step.}
#'   \item{`fn_full`}{A vector where each element is the objective function value at the given step.}
#' }
#' @export
adam <- function(par, fn, gr, nsteps,
                 alpha = 1e-3, beta = c(.9, .999), eps = 1e-8) {
  Fn <- rep(NA, nsteps)
  Par <- matrix(NA, nsteps, length(par))
  colnames(Par) <- names(par)
  mt <- 0
  vt <- 0
  thetat <- par
  betat <- c(1,1)
  for(tt in 1:nsteps) {
    gt <- gr(thetat)
    mt <- beta[1] * mt + (1-beta[1]) * gt
    vt <- beta[2] * vt + (1-beta[2]) * gt^2
    betat <- betat * beta
    mhat <- mt/(1-betat[1])
    vhat <- vt/(1-betat[2])
    thetat <- thetat - alpha * mhat / (sqrt(vhat) + eps)
    Par[tt,] <- thetat
    Fn[tt] <- fn(thetat)
  }
  imin <- which.min(Fn)
  list(par = Par[imin,], objective = Fn[imin],
       par_full = Par, Fn = Fn)
}
