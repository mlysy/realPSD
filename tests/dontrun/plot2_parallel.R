# to reproduce Figure 2 in the paper
require(realPSD)
# require(TMB)
require(tidyverse)
require(parallel)
# TODO: factor out an exportable function show_fit
# and perhaps a non-exported function show_fsim (f is for freq)
source("fitSHOW.R")
# data folder
data_path <- "~/Documents/data/R/realPSD/show_sim"

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
nsim <- 20
bin_size <- 100

# first, pregenerate exponentials
sim_expo <- FALSE
if(sim_expo) {
  for(ii in 1:nsim) {
    saveRDS(rexp(nf, rate = 1),
            file = file.path(data_path,
                             paste0("exp_sim_", ii, ".rds")))
  }
}

# now fitting
# this is for illustrative purposes mainly.
# i would probably do all the Q's and methods in the same loop,
# i.e., loop over only the datasets
fit_descr <- expand.grid(Q = Q_vec,
                         method = c("lp", "nls"),
                         data_id = 1:nsim,
                         stringsAsFactors = FALSE) %>% as_tibble()
nfit <- nrow(fit_descr)

success <- mclapply(1:nfit, function(ii) {
  # multi-assign elements of job: data_id, Q, method
  list2env(as.list(fit_descr[ii,]), envir = environment())
  # long form:
  ## data_id <- fit_descr$data_id[ii]
  ## Q <- fit_descr$Q[ii]
  ## method <- fit_descr$method[ii]
  psd <- readRDS(file.path(data_path,
                           paste0("exp_sim_", data_id, ".rds")))
  theta_hat <- tryCatch({
    fitSHOW(fseq, rfreq = psd,
            fs = fs, f0 = f0, Q = Q,
            k = k, T = Temp, Aw = Aw,
            binSize = bin_size, method = method)
  }, error = function(err) message("fitting error on job ", ii))
  if(!is.null(theta_hat)) {
    saveRDS(theta_hat,
            file = file.path(data_path,
                             paste0("show_fit_", ii, ".rds")))
    out <- TRUE
  } else out <- FALSE
  out
})

