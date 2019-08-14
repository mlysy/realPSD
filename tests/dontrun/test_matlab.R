# test optimization in MATLAB
require(realPSD)
require(tidyverse)
require(nlsr)
require(R.matlab)
source("fitSHOW_test.R")
# some initial constants
Time  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q_vec  <- c(1, 10, 100, 500)  # Quality factors
Q <- Q_vec[1]
k  <- 0.172               # Cantilever stiffness, N/m
Kb <- 1.381e-23           # Boltzmann's constant
Temp <- 298                  # Temperature, Kelvin
Aw <- 19000            # white noise, fm2/Hz 
Const <- 1e30    # 1m2 = 10^30 fm2
unit_conversion <- TRUE # if TRUE, convert m2/Hz to fm2/Hz
if(!unit_conversion) Aw <- Aw / Const          
# frequency domain
f_lb <- f0 - f0/sqrt(2) # frequency lower bound
f_ub <- f0 + f0/sqrt(2) # frequency upper bound
fseq <- seq(from = f_lb, to = f_ub, by = 1/Time) # frequency domain, Hz
nf <- length(fseq) # number of frequencies
# psd and periodogram values
psd <- psdSHO(fseq, f0, Q, k, Kb, Temp, unit_conversion) + Aw
set.seed(123)
sim_exp <- rexp(nf, rate = 1)
Y <- sim_exp * psd * fs
# save the periodogram ordinates (x) and frequency values (y) into MATLAB
# without binning beforehand as required by the MATLAB code
if(unit_conversion) {
  writeMat(con = file.path(paste0("./xPSD_Q", Q, "_unit_converted", ".mat")), 
        PSD_x = fseq)
  writeMat(con = file.path(paste0("./yPSD_Q", Q, "_unit_converted", ".mat")), 
        PSD_y = Y)
} else {
  writeMat(con = file.path(paste0("./xPSD_Q", Q, ".mat")), 
        PSD_x = fseq)
  writeMat(con = file.path(paste0("./yPSD_Q", Q, ".mat")), 
        PSD_y = Y)
}
# fit in R
nsim <- 1
fit_descr <- expand.grid(Q_level = Q,
                        method = c("nls", "lp", "mle"),
                        data_id = 1:nsim,
                        stringsAsFactors = FALSE)
# check the expanded grid
fit_descr
nfit <- nrow(fit_descr)
# fit all three methods
fit_data <- matrix(NA, nfit, 5)
bin_size <- 100
for(ii in 1:nfit) {
  # multi-assign elements of job: data_id, Q, method
  list2env(as.list(fit_descr[ii,]), envir = environment())
  # long form:
  data_id <- fit_descr$data_id[ii]
  Q <- fit_descr$Q_level[ii]
  method <- fit_descr$method[ii]
  fit_data[ii,] <- tryCatch({
    fitSHOW(fseq, sim_exp = sim_exp,
            f0 = f0, fs = fs, Q = Q,
            k = k, Temp = Temp, Aw = Aw,
            bin_size = bin_size, method = method,
            unit_conversion = unit_conversion)
  }, error = function(err) message("fitting error on job ", ii))
}
colnames(fit_data) <- c("As", "Q", "f0", "Aw", "k")
# display the result
fit_data <- fit_data %>% 
  as_tibble() %>% 
  add_column(method = c("nls", "lp", "mle")) %>%
  select(method, everything())
fit_data
