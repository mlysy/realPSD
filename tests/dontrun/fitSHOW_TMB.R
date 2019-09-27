fitSHOW_TMB <- function(fseq, Y, bin_size, method, phi, Temp, Kb) {
  # ---------- binning ----------
  fbar <- binning(fseq, bin_size = bin_size)
  Ybar <- binning(Y, bin_size = bin_size)
  Zbar <- log(Ybar)
  # bias_correction for LP method
  bias <- digamma(bin_size) - log(bin_size)
  # ---------- fitting ----------
  if (method == "lp") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_nlp",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gz <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "LP_zeta",
                                      fbar = matrix(fbar),
                                      Zbar = matrix(Zbar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
      zeta <- gz$fn(phi)
      # correct the bias
      exp(zeta - bias)
    }
  } else if (method == "nls") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_nlp",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "NLS_tau",
                                      fbar = matrix(fbar),
                                      Ybar = matrix(Ybar),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
      gt$fn(phi)
    }
  } else if (method == "mle") {
    obj <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_nlp",
                                      f = matrix(fseq),
                                      Y = matrix(Y),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
    get_tau <- function(phi) {
      gt <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                      method = "MLE_tau",
                                      f = matrix(fseq),
                                      Y = matrix(Y),
                                      fs = fs),
                          parameters = list(phi = matrix(0, 3, 1)),
                          silent = TRUE, DLL = "realPSD_TMBExports")
      gt$fn(phi)
    }
  } else {
    stop("method should be chosen from lp, nls and mle.")
  }
  # ---------- optimization -----------
  if(method == "mle" || method == "lp"){
    # ----- optim step-by-step -----
    start <- phi
    opt1 <- optim(start, fn = fn_fixed, 
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)],
      method = "Nelder-Mead") 
    start[2] <- opt1$par[2]
    opt2 <- optim(start, fn = fn_fixed,
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3],
      method = "Nelder-Mead")
    start[c(1,2)] <- opt2$par[c(1,2)]
    opt3 <- optim(start, fn = obj$fn)
    # opt3 <- optim(start, fn = obj$fn, gr = obj$gr, 
    #   method = "L-BFGS-B", lower = rep(0,3)) # may encounter error: L-BFGS-B needs finite values of 'fn'
    phi_hat <- opt3$par
    # start <- phi
    # opt1 <- pracma::fminsearch(f = fn_fixed, x0 = start, 
    #   obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)])
    # start[2] <- opt1$xmin[2]
    # opt2 <- pracma::fminsearch(f = fn_fixed, x0 = start, 
    #   obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3])
    # start[c(1,2)] <- opt2$xmin[c(1,2)]
    # opt3 <- pracma::fminsearch(f = obj$fn, x0 = start)
    # phi_hat <- opt3$xmin
    if(opt1$convergence != 0 || opt2$convergence != 0 || opt3$convergence != 0) phi_hat <- rep(NA, 3)
  } else {
    # ---------- lsqnonlin ---------
    # set some option parameters to avoid errors
    # tolx <- 1/(sum(phi^2)) # in order to compensate the squared norm of the parameter vector
    # tolg <- max(10*abs(obj$gr(phi)))
    tolx <- .Machine$double.eps^(2/3)
    tolg <- .Machine$double.eps^(2/3)
    # optimize Q (gamma), fix f0 and Rw
    opt1 <- pracma::lsqnonlin(fun = nls_res_fixed,
      x0 = phi, 
      options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)]) 
    # if(opt1$errno != 1) warning("NLS didn't converge at step 1.")
    # optimize Q and f0, fix Rw
    opt2 <- pracma::lsqnonlin(fun = nls_res_fixed, 
      x0 = opt1$x,
      options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3])
    # if(opt1$errno != 1) warning("NLS didn't converge at step 2.")
    # optimize all three parameters
    opt3 <- pracma::lsqnonlin(fun = nls_res, 
      x0 = opt2$x,
      options = list(tolx = tolx, tolg = tolg, maxeval = 1000),
      obj = obj)
    # opt3 <- optim(par = opt2$x, fn = obj$fn, gr = obj$gr, 
    #   method = "L-BFGS-B", lower = rep(0,3))
    # opt3 <- optim(par = opt2$x, fn = obj$fn, gr = obj$gr, 
    #   method = "BFGS")
    # if(opt1$errno != 1) warning("NLS didn't converge at step 3.")
    # return phi_hat
    # phi_hat <- opt3$par
    phi_hat <- opt3$x
  }
  # remove bad estimates 
  # if(any(phi_hat) <= 0) phi_hat <- rep(NA, 3)
  phi_hat[which(phi_hat <= 0)] <- NA
  # output
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2, unit should be the same as Aw
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2] # Q_hat
  param[3] <- Kb * Temp / (tau_hat * pi * param[1] * param[2]) # k_hat
  if(unit_conversion) {
    # param[4] <- phi_hat[3] * tau_hat * 1e30 # Aw_hat, same unit
    param[4] <- Kb * Temp / (pi * param[1] * param[2]) * 1e30
  } else {
    # param[4] <- phi_hat[3] * tau_hat
    param[4] <- Kb * Temp / (pi * param[1] * param[2])
  }
  names(param) <- c("f0_hat", "Q_hat", "k_hat", "Aw_hat")
  return(param)
}

# wrapper functions ----- method 1 ------
#' @param obj TMB obj
#' @param theta Parameter vector
#' @param fixed_flag Vector of TRUE/FALSE indicating which dimension of theta should be fixed
#' @param fixed_phi Vector of fixed values, length(fixed_phi) == length(which(fixed_id == TRUE))
fn_fixed <- function(theta, obj, fixed_flag, fixed_phi) {
  # set space without chaning the original theta
  theta_full <- rep(NA, length(theta))
  # fix part of theta
  # theta_full[which(fixed_flag == TRUE)] <- fixed_phi
  theta_full[fixed_flag == TRUE] <- fixed_phi
  # fill the remaining part with theta
  # theta_full[which(fixed_flag != TRUE)] <- theta
  theta_full[!fixed_flag] <- theta[!fixed_flag]
  # return
  obj$fn(theta_full)
}
gr_fixed <- function(theta, obj, fixed_flag, fixed_phi) {
  theta_full <- rep(NA, length(theta))
  theta_full[fixed_flag == TRUE] <- fixed_phi
  theta_full[!fixed_flag] <- theta[!fixed_flag]
  obj$gr(theta_full)
}
# wrapper function of the vector of residuals for NLS
nls_res <- function(phi, obj) {
  c(obj$simulate(phi)$RES)
} 
# wrapper of nls_res with fixed parameters
nls_res_fixed <- function(phi, obj, fixed_flag, fixed_phi) {
  phi_full <- rep(NA, length(phi))
  phi_full[fixed_flag == TRUE] <- fixed_phi
  phi_full[!fixed_flag] <- phi[!fixed_flag]
  nls_res(phi_full, obj)
}
