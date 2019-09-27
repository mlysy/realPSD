# generate plot3 
require(realPSD)
# ---------- SHO model parameters ----------
Time  <- 5                   # Total time, second
fs <- 1e7                    # Sampling frequency, Hz
f0 <- 33553                  # Resonance frequency, Hz
Q <- 100                     # Quality factors, unitless
k  <- 0.172                  # Cantilever stiffness, N/m
Kb <- 1.381e-23              # Boltzmann's constant, (m2*kg)/(s2*K)
Temp <- 298                  # Temperature, K
Const <- 1e30                # unit conversion, 1 m2 = 1e30 fm2
unit_conversion <- TRUE
Aw <- 19000                 # white noise, fm2/Hz 
if(!unit_conversion) Aw <- Aw / Const # if FALSE, we use the standard unit m2/Hz
# ---------- simulation ----------
fseq <- seq(from = 1/Time, to = fs - 1/Time, length.out = fs*Time) # frequency domain, Hz
N <- length(fseq)
cond <- which(fseq > f0-f0/sqrt(2) & fseq < f0+f0/sqrt(2))
fseq <- fseq[cond]
psd <- psdSHO(fseq, f0, Q, k, Temp, unit_conversion) + Aw
sin_fft <- fft_sin(fseq, N, f0, Q, fs, unit_conversion)
set.seed(2019)
n <- length(fseq)
x1 <- rnorm(n, 0, sqrt(1/2))
x2 <- rnorm(n, 0, sqrt(1/2))
sim_cnorm <- complex(real = x1, imaginary = x2)
Y <- sim_cnorm * sqrt(fs*psd)
Y <- (Y + sin_fft) * Conj(Y + sin_fft)
Y <- Re(Y)
pval <- 0.01
gstat <- Y / (fs * psd)
M <- max(gstat/sum(gstat))
fisherEq <- function(a, q, logSort, pval) {
  1 - fisherGstat(a, q, logSort) - pval
}
sol <- uniroot(fisherEq, c(0,M), 
  q = length(fseq), logSort = TRUE, pval = pval)
psd_line <- sol$root * psd * sum(gstat)
plot(x = fseq/10^3, y = Y/fs, 
  xlim = c((f0-1500)/10^3, (f0+1500)/10^3), 
  ylim = c(1e2, 1e10),
  xlab = "frequency (kHz)", ylab = expression(paste("PSD (", fm^2/Hz, ")")),
      type = "l", lwd = 0.3, log = "xy")
lines(x = fseq/10^3, y = psd_line, 
      col = "red")
legend("topright",
  legend = c("Baseline with electronic noise", "1% Denoising Threshold"), 
  col = c("black", "red"),
  bty = "n", # no box around the legend
  cex = 0.8, # size of the text,
  lwd = 1)
