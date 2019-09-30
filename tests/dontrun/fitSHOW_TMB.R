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
    exitflag <- 1 # set an exitflag: 1 means success, 0 means failure
    start <- phi # intial parameters
    # step 1: optimize Q
    opt1 <- optim(start, fn = fn_fixed, 
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)],
      method = "Nelder-Mead") 
    if(opt1$convergence != 0) exitflag <- 0
    # if(opt1$par[2] <= 0) { # if Nelder-Mead fails to find positive optim values
    #   opt11 <- optim(start, fn = fn_fixed, gr = gr_fixed,
    #     method = "BFGS",
    #     obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)])
    #   if(opt11$convergence != 0) exitflag <- 0
    #   start[2] <- opt11$par[2]
    # } else {
    #   start[2] <- opt1$par[2]
    # }
      start[2] <- opt1$par[2]
    # step 2: optimize f0, Q
    opt2 <- optim(start, fn = fn_fixed,
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3],
      method = "Nelder-Mead")
    if(opt2$convergence != 0) exitflag <- 0
    # if(any(opt2$par[c(1,2)] <= 0)) { # if Nelder-Mead fails to find positive optim values
    #   opt22 <- optim(start, fn = fn_fixed, gr = gr_fixed,
    #     method = "bfgs",
    #     obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3])
    #   if(opt22$convergence != 0) exitflag <- 0
    #   start[c(1,2)] <- opt22$par[c(1,2)]
    # } else {
    #   start[c(1,2)] <- opt2$par[c(1,2)]
    # }
      start[c(1,2)] <- opt2$par[c(1,2)]
    # step 3: optimize f0, Q, Rw all together
    opt3 <- optim(start, fn = obj$fn)
    if(opt3$convergence != 0) exitflag <- 0 
    # if(any(opt3$par <= 0)) {
    #   opt33 <- optim(start, fn = obj$fn, gr = obj$gr, method = "BFGS")
    #   if(opt33$convergence != 0) exitflag <- 0
    #   phi_hat <- opt33$par
    # } else{
    #   phi_hat <- opt3$par
    # }
      phi_hat <- opt3$par
  } else {
    # minpack.lm nls
    exitflag <- 1 # set an exitflag: 1 means success, 0 means failure
    start <- phi # initial param
    # step 1: optimize Q
    opt1 <- minpack.lm::nls.lm(
      par = start, 
      lower = rep(0,3),
      fn = nls_res_fixed, 
      obj = obj, fixed_flag = c(1,0,1), fixed_phi = phi[c(1,3)]) 
    # check if the optimization converged 
    # according to the docs of minpack.lm, 
    # info = 1,2,3,4 indicate a successful completion
    if(!is.element(opt1$info, c(1,2,3))) exitflag <- 0
    # update the initial param for the next step
    start[2] <- opt1$par[2]
    # step 2: optimize f0, Q
    opt2 <- minpack.lm::nls.lm(
      par = start, 
      lower = rep(0,3),
      fn = nls_res_fixed, 
      obj = obj, fixed_flag = c(0,0,1), fixed_phi = phi[3])  
    if(!is.element(opt2$info, c(1,2,3,4))) exitflag <- 0
    start[c(1,2)] <- opt2$par[c(1,2)]
    # step 3: optimize f0, Q, Rw
    opt3 <- minpack.lm::nls.lm(
      par = start,
      lower = rep(0,3),
      fn = nls_res,
      obj = obj)
    if(!is.element(opt3$info, c(1,2,3,4))) exitflag <- 0
    phi_hat <- opt3$par
  }
  # remove bad estimates 
  if(any(phi_hat <= 0)) phi_hat <- rep(NA, 3)
  if(exitflag != 1) phi_hat <- rep(NA, 3)
  # output
  tau_hat <- get_tau(phi_hat) # fitted tau = sigma^2, unit should be the same as Aw
  param <- rep(NA, 4) # allocate space for storage
  param[1] <- phi_hat[1] # f0_hat
  param[2] <- phi_hat[2] # Q_hat
  param[3] <- Kb * Temp / (tau_hat * pi * param[1] * param[2]) # k_hat
  if(unit_conversion) {
    param[4] <- phi_hat[3] * tau_hat * 1e30 # Aw_hat, same unit
  } else {
    param[4] <- phi_hat[3] * tau_hat
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
