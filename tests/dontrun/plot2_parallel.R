# to reproduce Figure 2 in the paper
require(realPSD)
# require(TMB)
require(parallel)
require(tidyverse)
require(tikzDevice)
# TODO: factor out an exportable function show_fit
# and perhaps a non-exported function show_fsim (f is for freq)
source("fitSHOW.R")
# data folder
# data_path_sim <- "~/Documents/data/R/realPSD/show_sim"
data_path_sim <- "~/realPSD/show_sim"
# data_path_fit <- "~/Documents/data/R/realPSD/show_fit"
data_path_fit <- "~/realPSD/show_fit"
# clear any existing files
# unlink(file.path(data_path_sim, "*"), recursive = TRUE) # we should keep the simulated expo rv's to save time
unlink(file.path(data_path_fit, "*"), recursive = TRUE)

# ---------- SHO model parameters ----------
Time  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q_vec  <- c(1, 10, 100, 500)  # Quality factors
k  <- 0.172               # Cantilever stiffness, N/m
# Kb <- 1.381e-23           # Boltzmann's constant
Temp <- 298                  # Temperature, Kelvin
Const <- 1e30
Aw <- 19000 / Const            # white noise, fm2/Hz to m2/Hz
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

# detect the number of cores
ncores <- detectCores()
# set seed for reproducibility
set.seed(2019, kind = "L'Ecuyer-CMRG")

# first, pregenerate exponentials
sim_expo <- TRUE
# system.time(
# if(sim_expo) {
#   for(ii in 1:nsim) {
#     saveRDS(rexp(nf, rate = 1),
#             file = file.path(data_path_sim,
#                             paste0("exp_sim_", ii, ".rds")))
#   }
# }
# )
# parallel version of generating random expo variables
message("\nTime spent on generating exponential random variables:\n")
system.time(
  if(sim_expo) {
    sim_success <- mclapply(1:nsim, function(ii) {
      tryCatch(
        saveRDS(rexp(nf, rate = 1),
            file = file.path(data_path_sim,
                            paste0("exp_sim_", ii, ".rds"))),
        error = function(err) 
          message(paste0("The ", ii, "-th simulation went wrong..."))
      )
      return(TRUE)
    }, mc.cores = ncores)
  }
)
# # check if all simulations were successful
# if(all(sim_success == TRUE)) {
#   message("Great! All simulations were successful!")
# } else {
#   err_index <- which(sim_success == FALSE)
#   message(paste0("The ", unname(err_index), "-th simulation(s) went wrong..."))
# }

# now fitting
# this is for illustrative purposes mainly.
# i would probably do all the Q's and methods in the same loop,
# i.e., loop over only the datasets
fit_descr <- expand.grid(Q_level = Q_vec,
                        method = c("nls", "lp", "mle"),
                        data_id = 1:nsim,
                        stringsAsFactors = FALSE)
nfit <- nrow(fit_descr)

# run the simulation
message("\nTime spent on fitting the parameters:\n")
system.time(
fit_success <- mclapply(1:nfit, function(ii) {
  # multi-assign elements of job: data_id, Q, method
  list2env(as.list(fit_descr[ii,]), envir = environment())
  # long form:
  data_id <- fit_descr$data_id[ii]
  Q <- fit_descr$Q_level[ii]
  method <- fit_descr$method[ii]
  r_exp <- readRDS(file.path(data_path_sim,
                          paste0("exp_sim_", data_id, ".rds")))
  theta_hat <- tryCatch({
    fitSHOW(fseq, sim_exp = r_exp,
            f0 = f0, fs = fs, Q = Q,
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
# # check if all iterations were successful
# if(all(fit_success == TRUE)) {
#   message("Great! All fitting jobs were successful!")
# } else {
#   err_index <- which(fit_success == FALSE)
#   message(paste0("The fitting job(s): ", unname(err_index), " had some errors..."))
# }
# all(fit_success == TRUE)

# ---------- Code below can be commented out for server computing ---------------
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
  mutate(method = factor(method, levels = c("nls", "lp", "mle"))) # factor column method 
# get a new dataset with ratios instead of fitted values
ratio_data <- fit_data %>% 
  mutate(f0_hat = f0_hat/f0) %>%
  mutate(k_hat = k_hat/k) %>%
  mutate(Aw_hat = Aw_hat/Aw) %>%
  mutate(Q_hat = case_when(
      Q_level == "Q = 1" ~ Q_hat/Q_vec[1],
      Q_level == "Q = 10" ~ Q_hat/Q_vec[2],
      Q_level == "Q = 100" ~ Q_hat/Q_vec[3],
      Q_level == "Q = 500" ~ Q_hat/Q_vec[4]
    )
  )
# boxplot
# Q_hat / Q
tikzDevice::tikz(file = "boxplot_Q.tex", width = 8, height = 2)
fig_Q <- ggplot(ratio_data, aes(x = Q_level, y = Q_hat, fill = method)) + 
  geom_boxplot(outlier.size = 0.8) +
  xlab(label = NULL)  +
  ylab(label = "$\\hat{Q}/Q$")
# ggsave("boxplot_Q.pdf")
print(fig_Q)
dev.off()

# k_hat / k
tikzDevice::tikz(file = "./boxplot_k.tex", width = 8, height = 2)
fig_k <- ggplot(ratio_data, aes(x = Q_level, y = k_hat, fill = method)) + 
  geom_boxplot(outlier.size = 0.8) +
  xlab(label = NULL)  +
  ylab(label = "$\\hat{k}/k$")
print(fig_k)
# ggsave("boxplot_k.pdf")
dev.off()

# f0_hat / f0
tikzDevice::tikz(file = "./boxplot_f0.tex", width = 8, height = 2)
fig_f0 <- ggplot(ratio_data, aes(x = Q_level, y = f0_hat, fill = method)) + 
  geom_boxplot(outlier.size = 0.8) + 
  xlab(label = NULL)  +
  ylab(label = "$\\hat{f_0}/f_0$")
print(fig_f0)
# ggsave("boxplot_f0.pdf")
dev.off()

# # arrange all the plots in one layout
# fig_all <- gridExtra::arrangeGrob(fig_Q, fig_k, fig_f0, ncol = 1)
# # ggsave(file = "boxplot_all.pdf", fig_all)

# ratio_Q500 <- ratio_data %>% filter(Q_level == "Q = 500")
# ggplot(ratio_Q500, aes(x = method, y = Q_hat)) + geom_boxplot(outlier.size = .8)


