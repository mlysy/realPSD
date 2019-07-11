# to reproduce Figure 2 in the paper
require(realPSD)
# ---------- SHO model parameters ----------
T  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q  <- c(1, 10, 100, 500)  # Quality factors
k  <- 0.172               # Cantilever stiffness, N/m
Kb <- 1.381e-23           # Boltzmann's constant
T <- 298                  # Temperature, Kelvin
Aw <- 19000               # white noise, fm2/Hz
sig2 <- Kb*T/(k*pi*f0*Q)  # variance sig2
Rw <- Aw/sig2
# alpha <- 0.55           # 1/f decay exponent
 
# ---------- simulate random datasets ----------
f_lb <- f0 - f0/sqrt(2) # frequency lower bound
f_ub <- f0 + f0/sqrt(2) # frequency upper bound
f <- seq(from = f_lb, to = f_ub, by = 1/T) # frequency domain, Hz
PSD <- matrix(NA, length(f), length(Q)) # each column corresponds to a Q factor
colnames(PSD) <- c("Q1", "Q10", "Q100", "Q500")
Y <- matrix(NA, length(f), length(Q))
colnames(Y) <- c("Q1", "Q10", "Q100", "Q500")
for(ii in 1:length(Q)) {
  PSD[,ii] <- psdSHO(f, f0, Q[ii], k, Kb, T, TRUE) + Aw
  Y[,ii] <- fs * PSD[, ii] * rexp(1, 1) # simulate the periodogram values
}
# ---------- binning ----------
binSize <- 100
nbins <- floor(length(f) / binSize)
# binning the frequency domain and periodogram values
fbar <- tapply(f, cut(f, nbins), mean)
Ybar <- apply(Y, 2, binning, nbins = nbins) 
Zbar <- log(Ybar)
# ---------- fitting ----------
tmod_lp <- TMB::MakeADFun(data = list(model_name = "SHOWFit",
                                       method = "LP_nlp",
                                       fbar = matrix(fbar),
                                       Zbar = matrix(Zbar)),
                           parameters = list(phi = matrix(c(f^2, f0*Q, Rw))),
                           silent = TRUE, DLL = "realPSD_TMBExports")
opt <- optim(tmod_lp$par, fn = tmod_lp$fn, gr = tmod_lp$gr, control = list(maxit = 1000))


