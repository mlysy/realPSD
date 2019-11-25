#--- develop and test stable fitting functions for each of the SHOW methods ----

library(realPSD)
library(TMB)
source("show-fitfunctions.R")

# SHOW model parameters
Time  <- 5                  # Total time, seconds
fs <- 1e7                   # Sampling frequency, Hz
f0 <- 33553                 # Resonance frequency, Hz
Q <- 500                    # Quality factor
k  <- 0.172                 # Cantilever stiffness, N/m
Temp <- 298                 # Temperature, K(elvin)
Sw <- 1.9e-26               # white noise power, m^2/Hz
Uconst <- 1e30              # convert m^2 to fm^2

# simulate dataset
# restricted frequency basis
frng <- f0 + c(-1,1) * f0/sqrt(2)
N <- Time*fs
fseq <- get_fseq(frng = frng, fs = fs, N = N)
nfreq <- length(fseq)
frng <- range(fseq)
# continuous-time psd (units of m^2/Hz)
cpsd <- show_psd(fseq = fseq, k = k, f0 = f0, Q = Q, Sw = Sw,
                 Temp = Temp)
# simulate PSD data directly in frequency domain
Xpsd <- sqrt(cpsd * fs/2) * (rnorm(nfreq) + 1i * rnorm(nfreq)) # complex normals
Ypsd <- abs(Xpsd)^2

# take a quick look
nplot <- 5000
iplot <- seq(1, nfreq, length = nplot)
plot(fseq[iplot]/1e3, Ypsd[iplot],
     xlab = "", ylab = "",
     pch = 16, cex = .5, log = "xy")
title(xlab = expression("Frequency "*(kHz)),
      ylab = expression("Power "*(m^2/Hz)),
      line = 2.5)
lines(fseq[iplot]/1e3, cpsd[iplot] * fs, lwd = 4, col = "red")

#--- LP estimator --------------------------------------------------------------

# let's start with the SHOW_log estimator

bin_size <- 100
fbar <- binning(fseq, bin_size = bin_size)
Zbar <- log(binning(Ypsd, bin_size = bin_size))
# normalize things for numerical stability
constZ <- mean(Zbar)

show_lp <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                       method = "LP_nlp",
                                       fbar = as.matrix(fbar),
                                       Zbar = as.matrix(Zbar - constZ),
                                       fs = fs/exp(bin_const(bin_size))),
                           parameters = list(phi = as.matrix(c(0,0,0))),
                           silent = TRUE, DLL = "realPSD_TMBExports")

# now equivalent to show_lp$simulate()$zeta (without AD)
## show_zeta <- TMB::MakeADFun(data = list(model = "SHOW_log",
##                                       method = "LP_zeta",
##                                       fbar = as.matrix(fbar),
##                                       Zbar = as.matrix(Zbar),
##                                       fs = fs/exp(bin_const(bin_size))),
##                           parameters = list(phi = as.matrix(c(0,0,0))),
##                           silent = TRUE, DLL = "realPSD_TMBExports")


# initial value
theta0 <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
phi0 <- log(theta0[1:3])


# direct quasi-newton method
lp_fit1 <- nlminb(start = phi0,
                  objective = show_lp$fn,
                  gradient = show_lp$gr)

# using optim
lp_fit2 <- optim(par = phi0,
                 fn = show_lp$fn,
                 gr = show_lp$gr,
                 method = "BFGS")
get_par(c(exp(lp_fit2$par),
          tau = exp(show_lp$simulate(lp_fit2$par)$zeta + constZ)),
        Temp = Temp)


# all in one
lp_fit2b <- show_fit_lp(f = fseq, Ypsd = Ypsd,
                        fs = fs, Temp = Temp,
                        bin_size = bin_size, phi0 = phi0,
                        fit_type = "incremental")
lp_fit2b$par

# using nlm
lp_fit3 <- nlm(f = function(phi) {
  out <- show_lp$fn(phi)
  attr(out, "gradient") <- show_lp$gr()
  out
}, p = phi0)

# using adam
lp_fit4 <- adam(theta0 = phi0,
                fn = show_lp$fn,
                gr = show_lp$gr,
                nsteps = 1e3)

plot(lp_fit4$Fn, log = "xy", type = "l")

# graphical check
library(optimCheck)

optim_proj(xsol = lp_fit1$par, fun = show_lp$fn)


show_psd_phi <- function(fseq, phi, Temp) {
  zeta <- show_lp$simulate(phi)$zeta
  par <- get_par(exp(c(phi, tau = zeta + constZ)), Temp = Temp)
  do.call(show_psd, args = c(list(fseq = fseq, Temp = Temp),
                             as.list(par)))
}

nplot <- 5000
iplot <- seq(1, nfreq, length = nplot)
plot(fseq[iplot],
     show_psd(fseq[iplot], k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp),
     xlab = "Frequency (Hz)", ylab = "Power (m^2/Hz)",
     type = "l", log = "xy")
lines(fseq[iplot],
      show_psd_phi(fseq[iplot], phi = lp_fit1$par, Temp = Temp),
      col = "red")


# single value optimization

# quick check that fn_fixed and gr_fixed work as expected
unlist(replicate(10, {
  fixed <- sample(c(TRUE, FALSE), 3, replace = TRUE)
  dfn <- fn_fixed(phi = phi0[!fixed], obj = show_lp,
                  fixed = fixed, phi0 = phi0) - show_lp$fn(phi0)
  dgr <- gr_fixed(phi = phi0[!fixed], obj = show_lp,
                  fixed = fixed, phi0 = phi0) - show_lp$gr(phi0)[!fixed]
  c(dfn, dgr)
}))

fixed <- c(FALSE, TRUE, TRUE)

optim(par = phi0[!fixed], fn = fn_fixed, gr = gr_fixed,
      method = "BFGS", obj = show_lp, fixed = fixed, phi0 = phi0)

#--- MLE estimator -------------------------------------------------------------

constY <- mean(Ypsd)

show_mle <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                       method = "MLE_nlp",
                                       f = as.matrix(fseq),
                                       Y = as.matrix(Ypsd/constY),
                                       fs = fs),
                           parameters = list(phi = as.matrix(c(0,0,0))),
                           silent = TRUE, DLL = "realPSD_TMBExports")

## show_tau <- TMB::MakeADFun(data = list(model = "SHOW_log",
##                                        method = "MLE_tau",
##                                        f = as.matrix(fseq),
##                                        Y = as.matrix(Ypsd),
##                                        fs = fs),
##                           parameters = list(phi = as.matrix(c(0,0,0))),
##                           silent = TRUE, DLL = "realPSD_TMBExports")

# initial value
theta0 <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
phi0 <- log(theta0[1:3])

# direct quasi-newton method
system.time({
  mle_fit1 <- nlminb(start = phi0,
                     objective = show_mle$fn,
                     gradient = show_mle$gr)
})

# using optim
system.time({
  mle_fit2 <- optim(par = phi0,
                    fn = show_mle$fn,
                    gr = show_mle$gr,
                    method = "BFGS")
})

get_par(c(exp(mle_fit2$par),
          tau = show_mle$simulate(mle_fit2$par)$tau * constY),
        Temp = Temp)


# all in one
system.time({
  mle_fit2b <- show_fit_mle(fseq = fseq, Ypsd = Ypsd,
                            fs = fs, Temp = Temp,
                            phi0 = phi0, fit_type = "incremental")
})
mle_fit2b$par

# using nlm
system.time({
  mle_fit3 <- nlm(f = function(phi) {
    out <- show_mle$fn(phi)
    attr(out, "gradient") <- show_mle$gr()
    out
  }, p = phi0)
})

oproj <- optim_proj(xsol = mle_fit1$par, fun = show_mle$fn)

plot(oproj, xnames = expression(log(f[0]), log(Q), log(R[w])))

show_psd_phi <- function(fseq, phi, Temp) {
  tau <- show_mle$simulate(phi)$tau
  par <- get_par(c(exp(phi), tau = tau * constY), Temp = Temp)
  do.call(show_psd, args = c(list(fseq = fseq, Temp = Temp),
                             as.list(par)))
}

nplot <- 5000
iplot <- seq(1, nfreq, length = nplot)
plot(fseq[iplot],
     show_psd(fseq[iplot], k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp),
     xlab = "Frequency (Hz)", ylab = "Power (m^2/Hz)",
     type = "l", log = "xy")
lines(fseq[iplot],
      show_psd_phi(fseq[iplot], phi = mle_fit2$par, Temp = Temp),
      col = "red")


#--- NLS estimator -------------------------------------------------------------

bin_size <- 100
fbar <- binning(fseq, bin_size = bin_size)
Ybar <- binning(Ypsd, bin_size = bin_size)
constY <-  mean(Ybar)

show_nls <- TMB::MakeADFun(data = list(model = "SHOW_log",
                                       method = "NLS_nlp",
                                       fbar = as.matrix(fbar),
                                       Ybar = as.matrix(Ybar/constY),
                                       fs = fs),
                           parameters = list(phi = as.matrix(c(0,0,0))),
                           silent = TRUE, DLL = "realPSD_TMBExports")

## show_tau <- TMB::MakeADFun(data = list(model = "SHOW_log",
##                                        method = "NLS_tau",
##                                        fbar = as.matrix(fbar),
##                                        Ybar = as.matrix(Ybar),
##                                        fs = fs),
##                           parameters = list(phi = as.matrix(c(0,0,0))),
##                           silent = TRUE, DLL = "realPSD_TMBExports")

# initial value
theta0 <- get_theta(k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp)
phi0 <- log(theta0[1:3])

# direct quasi-newton method
system.time({
  nls_fit1 <- nlminb(start = phi0,
                     objective = show_nls$fn,
                     gradient = show_nls$gr)
})

# using optim
system.time({
  nls_fit2 <- optim(par = phi0,
                    fn = show_nls$fn,
                    gr = show_nls$gr,
                    method = "BFGS")
})
get_par(c(exp(nls_fit2$par),
          tau = show_nls$simulate(nls_fit2$par)$tau * constY),
        Temp = Temp)


# all in one
system.time({
  fnscale <- show_nls$fn(phi0)
  nls_fit2b <- show_fit_nls(fseq = fseq, Ypsd = Ypsd,
                            fs = fs, Temp = Temp,
                            bin_size = bin_size, phi0 = phi0,
                            fit_type = "incremental",
                            control = list(fnscale = fnscale))
})
nls_fit2b$par


# using nlm
system.time({
  nls_fit3 <- nlm(f = function(phi) {
    out <- show_nls$fn(phi)
    attr(out, "gradient") <- show_nls$gr()
    out
  }, p = phi0)
})

# using adam
system.time({
nls_fit4 <- adam(theta0 = phi0,
                 fn = show_nls$fn,
                 gr = show_nls$gr,
                 nsteps = 1e3)
})

plot(nls_fit4$Fn)

oproj <- optim_proj(xsol = nls_fit2$par, fun = show_nls$fn, plot = FALSE)

plot(oproj, xnames = expression(log(f[0]), log(Q), log(R[w])))

show_psd_phi <- function(fseq, phi, Temp) {
  tau <- show_nls$simulate(phi)$tau
  par <- get_par(c(exp(phi), tau = tau * constY), Temp = Temp)
  do.call(show_psd, args = c(list(fseq = fseq, Temp = Temp),
                             as.list(par)))
}

nplot <- 5000
iplot <- seq(1, nfreq, length = nplot)
plot(fseq[iplot],
     show_psd(fseq[iplot], k = k, f0 = f0, Q = Q, Sw = Sw, Temp = Temp),
     xlab = "Frequency (Hz)", ylab = "Power (m^2/Hz)",
     type = "l", log = "xy")
lines(fseq[iplot],
      show_psd_phi(fseq[iplot], phi = nls_fit2$par, Temp = Temp),
      col = "red")
