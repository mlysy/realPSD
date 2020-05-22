#--- initialization of UFun object ---------------------------------------------

require(TMB)

model <- "UFun_Test"
compile(paste0(model, ".cpp"))
dyn.load(dynlib(model))

N <- 10
f <- matrix(rexp(N))

obj <- MakeADFun(data = list(f = f),
                 parameters = list(phi = matrix(rep(0, 3)),
                                   alpha = 0),
                 DLL = model, silent = TRUE)


phi <- matrix(rexp(3))
alpha <- rexp(1)
obj$simulate(c(phi, alpha))$U/obj$simulate(c(phi, 1))$U - alpha

## obj$simulate(alpha)$alpha - alpha
## range(obj$simulate(alpha)$f_out - f)
