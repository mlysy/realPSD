# ---------- simulation multiple OU processes and fit the model using TMB ----------
require(realPSD)
require(tidyverse)
require(TMB)
require(parallel)
source("fit-functions.R")
source("ou_sim.R")
# parallel computing settings
ncores <- 4
dir <- "./data"
dir_exp <- "./data_exp" # store simulated exponential rv
# unlink(file.path(dir, "*"), recursive = TRUE) # clean the data directory
# unlink(file.path(dir_exp, "*"), recursive = TRUE) # clean the data directory

# ---------- OU process parameter initialization ----------
alpha0 <- 2.3
beta0 <- sqrt(3.7)
theta0 <- c(alpha = alpha0, beta = beta0)
fmin = 1e-4
fmax = 1e2
fs <- 2 * fmax
Time = 1/fmin
dT <- 1/fs
N <- Time/dT
bin_size <- 100
fseq <- seq(from = 1/Time, to = fs/2, length.out = N/2)
# estimation freq range
frange <- c(1e-3,10)
findex <- c(which.min(abs(frange[1] - fseq)) : which.min(abs(frange[2] - fseq)))
nexp <- length(findex)
f <- fseq[findex]

# ---------- parallel simulation and fitting ----------
path <- 1:200
npath <- length(path)
descr_data <- expand.grid(method = c("NLS", "LP", "MLE"), path = path)
nfit <- nrow(descr_data)
descr_data <- cbind(job_id = 1:nfit, descr_data)

# parallel simulation
# system.time(
#   sim_success <- mclapply(1:npath, function(ii) {
#     # simulate an OU path
#     ou_obs <- ou_sim(gamma = alpha0, mu = 0, sigma = beta0, dt = dT, n_obs = N)
#     # calculate the periodogram
#     ou_psd <- periodogram(yTime = ou_obs, SF_s = fs, T_s = Time)
#     ou_psd <- ou_psd %>% as_tibble() %>% mutate(yFreq = yFreq/2) %>% rename(x = xFreq, y = yFreq)
#     # save the output and return a flag
#     ans <- tryCatch({
#       saveRDS(ou_obs, file = file.path(dir, paste0("ou_obs_path=", ii, ".rds")))
#       saveRDS(ou_psd, file = file.path(dir, paste0("ou_psd_path=", ii, ".rds")))
#       TRUE
#     }, error = function(err) FALSE
#     )
#     return(ans)
#   }, mc.cores = ncores)
# )
# parallel simulation of exponential random variables (Prop 1)
sim_expo <- TRUE
message("\nTime spent on generating exponential random variables:\n")
system.time(
  if(sim_expo) {
    sim_success <- mclapply(1:npath, function(ii) {
      tryCatch(
        saveRDS(rexp(nexp, rate = 1),
            file = file.path(dir_exp,
                            paste0("exp_sim_path=", ii, ".rds"))),
        error = function(err)
          message(paste0("The ", ii, "-th simulation went wrong..."))
      )
      return(TRUE)
    }, mc.cores = ncores)
  }
)
# check sim_success
which(sim_success == FALSE)

# parallel fitting
system.time(
  fit_success <- mclapply(1:nfit, function(ii) {
    # assign some const variables
    list2env(as.list(descr_data[ii,]), envir = environment())
    method <- descr_data$method[ii]
    path_id <- descr_data$path[ii]
    # load the simulated data
    # ou_psd <- readRDS(file.path(dir, paste0("ou_psd_path=", path_id, ".rds")))
    # xPSD <- ou_psd$x
    # yPSD <- ou_psd$y
    # set the freq range before fitting
    # frange <- c(1e-3,10)
    # findex <- c(which.min(abs(frange[1] - xPSD)) : which.min(abs(frange[2] - xPSD)))
    # f <- xPSD[findex]
    # Y <- yPSD[findex]

    # load sim exp
    r_exp <- readRDS(file.path(dir_exp,
                          paste0("exp_sim_path=", path_id, ".rds")))
    Y <- r_exp * beta0^2 / (4*pi^2*f^2 + alpha0^2)
    # plot(f, Y, type = "l", log = "xy")
    # curve(beta0^2 * ou_psd_r(f = x, phi = log(theta0)),
    #   col = "red", add = TRUE)
    # fitting the OU model using TMB
    model <- "ou_log"
    suppressWarnings({
      # compile(paste0(model, ".cpp"), PKG_CXXFLAGS = paste0("-I", system.file("include", package = "realPSD")))
      dyn.load(dynlib(model))
    })
    # initial value
    phi0 <- log(c(alpha = 1))
    # NLS
    if(method == "NLS") {
      fit_nls <- ou_fit_nls(model = model, fseq = f, Ypsd = Y, fs = fs,
                      bin_size = bin_size, bin_type = "mean",
                      phi0 = phi0, optimizer = "optim",
                      vcov = TRUE, get_jac = TRUE)
      # est
      param <- fit_nls$par
      # se
      # N_bin <- floor(length(Y)/bin_size)
      nls_jac <- fit_nls$jac
      B <- t(nls_jac) %*% nls_jac
      se <- setNames(sqrt(diag(fit_nls$cov %*% B %*% fit_nls$cov)), nm = c("alpha", "beta"))
    } else if(method == "LP") {
      fit_lp <- ou_fit_lp(model = model, fseq = f, Ypsd = Y, fs = fs,
                      bin_size = bin_size, bin_type = "mean",
                      phi0 = phi0, optimizer = "optim",
                      vcov = TRUE)
      param <- fit_lp$par
      se <- setNames(sqrt(diag(2/bin_size * fit_lp$cov)), nm = c("alpha", "beta"))
    } else if(method == "MLE") {
      fit_mle <- ou_fit_mle(model = model, fseq = f, Ypsd = Y, fs = fs,
                      bin_size = bin_size, bin_type = "mean",
                      phi0 = phi0, optimizer = "optim",
                      vcov = TRUE)
      param <- fit_mle$par
      se <- setNames(sqrt(diag(fit_mle$cov)), nm = c("alpha", "beta"))
    } else message(paste0("Unknown method for job = ", ii))
    # save the fitted results and return a flag
    ans <- tryCatch(
      {
        saveRDS(param, file = file.path(dir, paste0("param-job_id=", ii, ".rds")))
        saveRDS(se, file = file.path(dir, paste0("se-job_id=", ii, ".rds")))
        TRUE
      }, error = function(err) FALSE
    )
    return(ans)
  }, mc.cores = ncores)
)
# check if there is any failure
which(fit_success == FALSE)

# ---------- summary ----------
est_list <- list()
se_list <- list()
for(ii in 1:nfit) {
  est_list[[ii]] <- readRDS(file.path(dir, paste0("param-job_id=", ii, ".rds")))
  se_list[[ii]] <- readRDS(file.path(dir, paste0("se-job_id=", ii, ".rds")))
}
est_data <- do.call(rbind, est_list)
colnames(est_data) <- c("alpha_est", "beta_est")
summary <- cbind(descr_data[,c("method", "path")], est_data)
se_data <- do.call(rbind, se_list)
colnames(se_data) <- c("alpha_se", "beta_se")
summary <- cbind(summary, se_data)
summary <- summary %>% as_tibble() %>% relocate(alpha_se, .after = alpha_est)
summary %>% print(n=100)
# verify mean(se^2) == var(est)
var_tbl <- summary %>% group_by(method) %>%
  summarize(
    # alpha_var = mean((alpha_est-alpha0)^2),
    alpha_var = var(alpha_est),
    alpha_se2 = mean(alpha_se^2),
    # beta_var = mean((beta_est-beta0)^2),
    beta_var = var(beta_est),
    beta_se2 = mean(beta_se^2)
  ) %>%
  mutate(rel_diff_alpha = (alpha_se2 - alpha_var)/alpha_var) %>%
  mutate(rel_diff_beta = (beta_se2 - beta_var)/beta_var) %>%
  relocate(rel_diff_alpha, .after = alpha_se2) %>% print(n=100)


#--- scratch -------------------------------------------------------------------

obj_fun <- function(theta, Ybar, fbar) {
  alpha <- theta[1]
  beta <- theta[2]
  psd <- beta^2 / (4*pi^2*fbar^2 + alpha^2)
  (Ybar - psd)^2
}

theta_hat <- fit_nls$par
bin_size <- 100
Ybar <- binning(Y, bin_size, "mean")
fbar <- binning(f, bin_size, "mean")
constY <- mean(Ybar)

he <- numDeriv::hessian(func = function(theta) sum(obj_fun(theta, Ybar, fbar)),
                        x = theta_hat)
vcov <- chol2inv(chol(he))

jac <- numDeriv::jacobian(func = function(theta) obj_fun(theta, Ybar, fbar),
                          x = theta_hat)

vcov %*% crossprod(jac) %*% vcov
fit_nls$cov %*% crossprod(fit_nls$jac) %*% fit_nls$cov

head(obj_fun(fit_nls$par, Ybar, fbar))
