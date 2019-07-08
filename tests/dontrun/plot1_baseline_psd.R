# plot Figure 1 in the paper(latest version)
# total time
T_s <- 5 
# sampling frequency
SF_s <- 1e6 
# resonance frequency in Hz
f0_s <- 33553 
# quality factors
Q1 <- 1 
Q10 <- 10 
Q100 <- 100 
Q500 <- 500 
# 1/f decay exponent
alpha_s <- 0.55 
# cantilever stiffness N/m
k_s <- 0.172 
# Boltzmann's constant
Kb <- 1.381e-23 
# Kelvin
T <- 298 
# unit conversion
CONST <- 1e30   
# white noise
Aw_s <- 19000  
# frequency domain (frequency grid)
f <- seq(from = 1/T_s, to = SF_s, length.out = SF_s*T_s)
# x-axis for the PSD plot
xPSD <- f
# PSD of the different SHO models with Q1, Q10, Q100, Q500
Q <- c(Q1, Q10, Q100, Q500)
yPSD <- matrix(NA, length(f), length(Q)) # each column corresponds to a Q factor
# when we assign the values, we also convert the unit of yPSD
for(ii in 1:length(Q)) {
  yPSD[,ii] <- psd_sho(f, f0_s, Q[ii], k_s, Kb, T)*CONST + Aw_s
}
# change x-axis unit to kHz
xPSD <- xPSD / 1000
# plot Figure 2
# set a palette with black(#000000)
cbp <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#CC79A7", "#000000")
matplot(xPSD, yPSD, log = "xy", xlim = c(10, 55), ylim = c(1e4, 1e9),
  type = "l", lty = 1, col = cbp[4:1])
# log10.axis(side = 2, at = seq(4,9,1))
# minor.ticks.axis(2, 9, mn = 4, mx = 9)
axis(side = 1, at = seq(10, 55, 5))
# marks <- c(1e4, 1e5, 1e6, 1e7, 1e8, 1e9)
# format(y,scientific=FALSE)
# axis(2,at=marks,labels=format(marks,scientific=FALSE))
