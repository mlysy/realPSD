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
# yes, unit conversion
CONST <- 1e30   
# white noise
Aw_s <- 19000  
# 1/f noise
Af_s <- 0 
# frequency domain
seq(from = 1/T_s, to = SF_s, length.out = SF_s*T_s); % freq space.


# power spectral density of SHO model
