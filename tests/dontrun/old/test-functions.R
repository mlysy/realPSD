# test helper functions
# simulation functions
sim_f <- function(n) runif(n, 0, 2*n)
sim_phi <- function() c(f0 = runif(1, 1e3, 1e5), gamma = rexp(1)*1e5, Rw = rexp(1)/1e5)
