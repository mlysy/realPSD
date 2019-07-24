# to reproduce Figure 2 in the paper
require(realPSD)
# require(TMB)
# require(tidyverse)
require(parallel)
# TODO: factor out an exportable function show_fit
# and perhaps a non-exported function show_fsim (f is for freq)
source("fitSHOW.R")
# data folder
# data_path_sim <- "~/Documents/data/R/realPSD/show_sim"
data_path_sim <- "~/realPSD/show_sim"
# data_path_fit <- "~/Documents/data/R/realPSD/show_fit"
data_path_fit <- "~/realPSD/show_fit"
# clear any existing files
unlink(file.path(data_path_sim, "*"), recursive = TRUE)
unlink(file.path(data_path_fit, "*"), recursive = TRUE)

# ---------- SHO model parameters ----------
Time  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q_vec  <- c(1, 10, 100, 500)  # Quality factors
k  <- 0.172               # Cantilever stiffness, N/m
# Kb <- 1.381e-23           # Boltzmann's constant
Temp <- 298                  # Temperature, Kelvin
Aw <- 19000               # white noise, fm2/Hz
# sig2 <- Kb*T/(k*pi*f0*Q)  # variance sig2
# Rw <- Aw/sig2
# alpha <- 0.55           # 1/f decay exponent

# ---------- simulate random datasets ----------
f_lb <- f0 - f0/sqrt(2) # frequency lower bound
f_ub <- f0 + f0/sqrt(2) # frequency upper bound
fseq <- seq(from = f_lb, to = f_ub, by = 1/Time) # frequency domain, Hz
nf <- length(fseq) # number of frequencies

# ---------- simulation ----------
nsim <- 1000
bin_size <- 100

# first, pregenerate exponentials
sim_expo <- TRUE
if(sim_expo) {
  for(ii in 1:nsim) {
    saveRDS(rexp(nf, rate = 1),
            file = file.path(data_path_sim,
                            paste0("exp_sim_", ii, ".rds")))
  }
}

# now fitting
# this is for illustrative purposes mainly.
# i would probably do all the Q's and methods in the same loop,
# i.e., loop over only the datasets
fit_descr <- expand.grid(Q_level = Q_vec,
                        method = c("lp", "nls", "mle"),
                        data_id = 1:nsim,
                        stringsAsFactors = FALSE)
nfit <- nrow(fit_descr)

# detect the number of cores
ncores <- detectCores()
# set seed for reproducibility
set.seed(123, kind = "L'Ecuyer-CMRG")
# run the simulation
system.time(
success <- mclapply(1:nfit, function(ii) {
  # multi-assign elements of job: data_id, Q, method
  list2env(as.list(fit_descr[ii,]), envir = environment())
  # long form:
  data_id <- fit_descr$data_id[ii]
  Q <- fit_descr$Q[ii]
  method <- fit_descr$method[ii]
  r_exp <- readRDS(file.path(data_path_sim,
                          paste0("exp_sim_", data_id, ".rds")))
  theta_hat <- tryCatch({
    fitSHOW(fseq, sim_exp = r_exp,
            fs = fs, f0 = f0, Q = Q,
            k = k, Temp = Temp, Aw = Aw,
            bin_size = bin_size, method = method)
  }, error = function(err) message("fitting error on job ", ii))
  if(!is.null(theta_hat)) {
    saveRDS(theta_hat,
            file = file.path(data_path_fit,
                            paste0("show_fit_", ii, ".rds")))
    out <- TRUE
  } else out <- FALSE
  out
}, mc.cores = ncores)
)
# check if all iterations were successful
all(success == TRUE)

# ---------- Code below should be commented out for server computing ---------------
# read simulated data into workspace
for(ii in 1:nfit) {
  assign(paste0("fit_", ii), 
        readRDS(file.path(data_path_fit, paste0("show_fit_", ii, ".rds"))))
}
# combine each iteration vector into a data frame
fit_list <- list()
for(ii in 1:nfit) {
  fit_list[[ii]] <- readRDS(file.path(data_path_fit, paste0("show_fit_", ii, ".rds"))) 
}
fit_data <- do.call(rbind, fit_list)
fit_data <- cbind(fit_data, fit_descr[,c("Q_level", "method")])
# then manipulate the data frame using tidyverse toolbox
fit_data <- fit_data %>% as_tibble() %>%
  mutate(Q_level = factor(Q_level,  # convert the column Q_level into a factor
    levels = c(1,10,100,500), labels = c("Q = 1", "Q = 10", "Q = 100", "Q = 500"))) %>% 
  mutate(method = factor(method, ordered = FALSE)) # factor column method 
# get a new dataset with ratios instead of fitted values
ratio_data <- fit_data %>% 
  mutate(f0_hat = f0_hat/f0) %>%
  mutate(k_hat = k_hat/k) %>%
  mutate(Aw_hat = Aw_hat/Aw) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 1", Q_hat/Q_vec[1], Q_hat)) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 10", Q_hat/Q_vec[2], Q_hat)) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 100", Q_hat/Q_vec[3], Q_hat)) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 500", Q_hat/Q_vec[4], Q_hat))
  # boxplot
# Q_hat / Q
ggplot(ratio_data, aes(x = Q_level, y = Q_hat, fill = method)) + geom_boxplot()
# k_hat / k
ggplot(ratio_data, aes(x = Q_level, y = k_hat, fill = method)) + geom_boxplot()
# f0_hat / f0
ggplot(ratio_data, aes(x = Q_level, y = f0_hat, fill = method)) + geom_boxplot()
