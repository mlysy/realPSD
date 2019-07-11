# to reproduce Figure 2 in the paper
# ---------- SHO model parameters ----------
T  <- 5                   # Total time, second
fs <- 1e7                 # Sampling frequency, Hz
f0 <- 33553               # Resonance frequency, Hz
Q  <- c(1, 10, 100, 500)  # Quality factor
k  <- 0.172               # Cantilever stiffness, N/m
Kb <- 1.381e-23           # Boltzmann's constant
T <- 298                  # Temperature, Kelvin
# alpha <- 0.55           # 1/f decay exponent
 
# ---------- simulate random datasets ----------
f <- seq(from = 1/T, to = fs, length.out = fs*T) # frequency domain


# ---------- binning ----------
