require(realPSD)
require(TMB)
require(parallel)
source("show-fitfunctions.R")

# SHOW model parameters
Time  <- 5                  # Total time, seconds
fs <- 1e7                   # Sampling frequency, Hz
f0 <- 33553                 # Resonance frequency, Hz
Q_vals <- c(1,10,100,500)   # Quality factors (unitless)
k  <- 0.172                 # Cantilever stiffness, N/m
Temp <- 298                 # Temperature, K(elvin)
Sw <- 1.9e-26               # white noise power, m^2/Hz

# other parameters
output_path <- "show-fit_fig2_sim1-" # path and file prefix to save output
fnrg <- f0 + c(-1, 1) * f0/sqrt(2) # restricted frequency basis
nsim <- 100                 # number of simulated datasets
bin_size <- 100 # bin size
ncores <- detectCores(logical = FALSE) # number of cores to use

job_descr <- expand.grid(Q = Q_vals,
                         data_id = 1:nsim,
                         stringsAsFactors = FALSE)
njobs <- nrow(job_descr)
rownames(job_descr) <- 1:njobs

# pre-generate random seeds to make things fully reproducible
RNGkind("L'Ecuyer-CMRG")
set.seed(2020)
# assigna a random seed to each dataset
Random_Seeds <- matrix(NA, nrow = nsim, ncol = length(.Random.seed))
s <- .Random.seed
for(ii in 1:nsim) {
  Random_Seeds[ii,] <- s
  s <- nextRNGStream(s)
}

# function to execute each job
job_fun <- function(job_id, save = FALSE) {
  # extract job parameters
  Q <- job_descr$Q[job_id]
  data_id <- job_descr$data_id[job_id]
  #----- simulate data -----
  .Random.seed <<- Random_Seeds[data_id,]
  # frequency basis
  N <- Time * fs
  fseq <- get_fseq(frng = frng, fs = fs, N = N)
  nfreq <- length(fseq)
  frng <- range(fseq)
  # continuous-time psd (units of m^2/Hz)
  cpsd <- show_psd(fseq = fseq, k = k, f0 = f0, Q = Q, Sw = Sw,
                   Temp = Temp)
  # simulate PSD data directly in frequency domain
  Xpsd <- sqrt(cpsd * fs/2) * (rnorm(nfreq) + 1i * rnorm(nfreq))
  Ypsd <- abs(Xpsd)^2
  # ----- fitting -----
  # initial value
  theta0 <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
  phi0 <- log(theta0[1:3])
  # lp method
  lp_fit <- tryCatch({
    show_fit_lp(f = fseq, Ypsd = Ypsd,
                        fs = fs, Temp = Temp,
                        bin_size = bin_size, phi0 = phi0,
                fit_type = "direct")
  }, error = catch_error)
  # nls method
  nls_fit <- tryCatch({
    show_fit_nls(f = fseq, Ypsd = Ypsd,
                 fs = fs, Temp = Temp,
                 bin_size = bin_size, phi0 = phi0,
                 fit_type = "direct")
  }, error = catch_error)
  # mle method
  mle_fit <- tryCatch({
    show_fit_mle(f = fseq, Ypsd = Ypsd,
                 fs = fs, Temp = Temp,
                 bin_size = bin_size, phi0 = phi0,
                 fit_type = "direct")
  }, error = catch_error)
  out <- list(lp = lp_fit, nls = nls_fit, mle = mle_fit)
  if(save) {
    saveRDS(out, file = paste0(output_path, job_id, ".rds"))
  } else {
    return(out)
  }
}

## # test
## out <- job_fun(3)
## out2 <- readRDS("show-fit_fig2_sim1-3.rds")

# set up cluster
cl <- makeCluster(spec = ncores)
clusterExport(cl,
              varlist = c("Time", "fs", "f0", "Q_vals", "k", "Temp", "Sw",
                          "output_path", "frng", "bin_size",
                          "job_descr", "Random_Seeds"))
clusterEvalQ(cl, expr = {
  require(realPSD)
  require(TMB)
  source("show-fitfunctions.R")
})

system.time({
  parLapply(cl, X = 1:4, fun = job_fun, save = TRUE)
})

stopCluster(cl) # shut down cores
