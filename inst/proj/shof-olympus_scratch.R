#' # Scratch
#'
#' This code won't be run, nor can it be expected to work as active code is updated.  It is kept here solely for reference.
#'
#' ## Plot Data

#+ plot_data, eval = FALSE
## plot_psd <- function(psd_emp, bin_size = 100, bin_type = "mean",
##                      log = "xy") {
##   psd_bin <- tibble(fbar = binning(psd_emp$freq,
##                                    bin_size = bin_size,
##                                    bin_type = bin_type),
##                     Ybar = binning(psd_emp$Ypsd,
##                                    bin_size = bin_size,
##                                    bin_type = bin_type))
##   plt <- psd_bin %>%
##     ggplot(aes(x = fbar, y = Ybar)) +
##     geom_line(size = 0.3, color = "deepskyblue1") +
##     xlab(expression(Frequency (kHz))) +
##     ylab(expression(PSD (m^2/Hz))) +
##     theme(axis.text = element_text(size=10),
##           axis.title = element_text(size=10,face="plain"))
##   if(log %in% c("x", "xy")) {
##     plt <- plt + scale_x_log10()
##   }
##   if(log %in% c("y", "xy")) {
##     plt <- plt + scale_y_log10()
##   }
##   ## geom_line() +
##   ##   scale_x_log10() +
##   ##   scale_y_log10()
## }

## psd_emp %>%
## clrs <- c(brewer.pal(3, "Set1")[2], "black")
psd_color <- "deepskyblue1"
tick_breaks <- 10^(-30:30)
minor_breaks <- as.matrix(expand.grid(1:9, tick_breaks))
minor_breaks <- minor_breaks[,1] * minor_breaks[,2]
label10 <- function(x) parse(text = paste0("10^", log10(x)))
line_size <- .3
psd_bin <- sapply(psd_emp, binning,
                  bin_size = 100, bin_type = "mean") %>%
  as_tibble()
## psd_bin <- readRDS("psd_bin.rds")
plt_main <- psd_bin %>%
  transmute(x = freq/1000, y = Ypsd) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(size = line_size, color = psd_color) +
  scale_x_log10(# trans = log10_trans(),
    breaks = tick_breaks,
    minor_breaks = minor_breaks,
    labels = label10
  ) +
  scale_y_log10(
    ## breaks = tick_breaks
    minor_breaks = 10^(-30:30),
    labels = label10
  ) +
  ## xlab(expression("Frequency (Hz)")) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("(a) Periodogram") +
  ## ylab(expression("PSD "*(m^2/Hz))) +
  annotation_logticks(size = 0.3) +
  theme_bw()
  ## theme(
  ##   plot.title = element_text(
  ##     ## margin = margin(t = 20, b = -30),
  ##     ## vjust = -1,
  ##     ## vjust = -10,
  ##     hjust = .5
  ##   )
  ## )
  ## theme_min2()
  ## theme(axis.text = element_text(size=10),
  ##       axis.title = element_text(size=10,face="plain"))
## plt_main

## frng1 <- c(20, 60)*1000
## frng2 <- c(100, 300)*1000
## plt_eig1 <- psd_emp %>%
##   filter(frng1[1] <= freq & freq <= frng1[2]) %>%
##   sapply(binning,
##          bin_size = 100, bin_type = "mean") %>%
##   as_tibble() %>%
##   transmute(x = freq/1000, y = Ypsd)
## plt_eig2 <- psd_emp %>%
##   filter(frng2[1] <= freq & freq <= frng2[2]) %>%
##   sapply(binning,
##          bin_size = 100, bin_type = "mean") %>%
##   as_tibble() %>%
##   transmute(x = freq/1000, y = Ypsd, mode = 2)
## plt_eigs <- bind_rows(plt_eig1, plt_eig2)


frange <- c(20, 60)*1000
plt_eig1 <- psd_bin %>%
  filter(frange[1] <= freq & freq <= frange[2]) %>%
  transmute(x = freq/1000, y = Ypsd, mode = 1)
plt_eig1 <- plt_eig1 %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(size = line_size, color = psd_color) +
  ## scale_x_log10(breaks = base_breaks(5)) +
  scale_y_log10(
    minor_breaks = minor_breaks,
    labels = label10
  ) +
  ## xlab(expression("Frequency (KHz)")) +
  ## ylab(expression("PSD "*(m^2/Hz))) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("(b) 1st Eigenmode") +
  annotation_logticks(size = 0.3, sides = "l") +
  theme_bw()
  ## theme(axis.text = element_text(size=10),
  ##       axis.title = element_text(size=10,face="plain"))
## plt_eig1

frange <- c(100, 300)*1000
plt_eig2 <- psd_bin %>%
  filter(frange[1] <= freq & freq <= frange[2]) %>%
  transmute(x = freq/1000, y = Ypsd)
plt_eig2 <- plt_eig2 %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(size = line_size, color = psd_color) +
  ## scale_x_log10(breaks = base_breaks(5)) +
  scale_y_log10(
    minor_breaks = minor_breaks,
    labels = label10
  ) +
  ## xlab(expression("Frequency (KHz)")) +
  ## ylab(expression("PSD "*(m^2/Hz))) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("(c) 2nd Eigenmode") +
  annotation_logticks(size = 0.3, sides = "l") +
  theme_bw()
  ## theme(axis.text = element_text(size=10),
  ##       axis.title = element_text(size=10,face="plain"))
## plt_eig2

plot_group <- grid.arrange(
  grobs = list(plt_main, plt_eig1, plt_eig2),
  layout_matrix = rbind(c(1,1), c(2,3)),
  left = grid::textGrob(expression("PSD "*(m^2/Hz)),
                        rot = 90),
  bottom = grid::textGrob(expression("Frequency (KHz)"))
)
plot_group

## plot_group <- ggarrange(
##   plt_main,
##   ggarrange(plt_eig1, plt_eig2, ncol = 2,
##             labels = c("(b) 1st Eigenmode", "(c) 2nd Eigenmode"),
##             hjust = -4, vjust = 4,
##             font.label = list(size = 10, face = "plain")),
##   nrow = 2,
##   labels = c("(a) Periodogram", ""),
##   hjust = -11, vjust = 4,
##   font.label = list(size = 10, face = "plain")
## ) +
##   ylab("PSD")
## plot_group

#' ## Fit NLS Estimator

#+ nls_fit_test, eval = FALSE
# frequency range
frange <- c(20, 60)*1000
bin_size <- 100
psd_data <- psd_emp %>% filter(freq >= frange[1] & freq <= frange[2])

# preliminary fit
set.seed(97)
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
shof_obj$set_est(est_type = "nls",
                 bin_size = bin_size, bin_type = "median")

nlp <- shof_obj$nlp()
fit <- optim(par = phi0, fn = nlp$fn, gr = nlp$gr, method = "BFGS")

# check fit
par(mfrow = c(1,1))
plot(shof_obj$fbar, shof_obj$Ybar / log(2), type = "l", log = "xy")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")

# denoising
theta_fit <- shof_obj$to_theta(phi = fit$par,
                               zeta = shof_obj$zeta(fit$par)$fn())
psd <- shof_psd(freq = psd_data$freq, theta = theta_fit, Temp = Temp)
Ypsd_den <- psd_denoise(Ypsd = psd_data$Ypsd, alpha = .01,
                        psd = psd, method = "fisherG")
Ypsd_den2 <- psd_denoise(Ypsd = psd_data$Ypsd, alpha = .01,
                         psd = psd, method = "BH")

bin_type <- "mean"
fbar <- binning(psd_data$freq, bin_size, bin_type = bin_type)
Ybar <- binning(psd_data$Ypsd, bin_size, bin_type = bin_type)
Ybar_den <- binning(Ypsd_den, bin_size, bin_type = bin_type)
Ybar_den2 <- binning(Ypsd_den2, bin_size, bin_type = bin_type)

par(mfrow = c(1,3))
xlim <- frange
type <- "l"
plot(fbar, Ybar, type = type, log = "xy",
     xlim = xlim, main = "No Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")
plot(fbar, Ybar_den, type = type, log = "xy",
     xlim = xlim, main = "fisherG Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")
plot(fbar, Ybar_den2, type = type, log = "xy",
     xlim = xlim, main = "BH Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")

# final fit
shof_obj$set_psd(freq = psd_data$freq, Ypsd = Ypsd_den)
shof_obj$set_est(est_type = "nls",
                 bin_size = bin_size, bin_type = "mean")

nlp <- shof_obj$nlp()
fit <- optim(par = fit$par, fn = nlp$fn, gr = nlp$gr, method = "BFGS")
nls_est <- to_est(fit$par, shof_obj)

# short version
set.seed(97)
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
nls_fit <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                   est_type = "nls", bin_size = bin_size)

# check
range(nls_fit$coef - nls_est$coef)
range(nls_fit$se - nls_est$se)
range((nls_fit$psd$Ypsd_fg - Ypsd_den)/Ypsd_den)


#' ## Fit LP Estimator

#+ lp_fit_test, eval = FALSE
# frequency range
frange <- c(20, 60)*1000
bin_size <- 100
psd_data <- psd_emp %>% filter(freq >= frange[1] & freq <= frange[2])

# preliminary fit
set.seed(300)
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
shof_obj$set_est(est_type = "lp",
                 bin_size = bin_size, bin_type = "median")
nlp <- shof_obj$nlp()
fit <- optim(par = phi0, fn = nlp$fn, gr = nlp$gr, method = "BFGS")

# denoising
theta_fit <- shof_obj$to_theta(phi = fit$par,
                               zeta = shof_obj$zeta(fit$par)$fn())
psd <- shof_psd(freq = psd_data$freq, theta = theta_fit, Temp = Temp)
Ypsd_den <- psd_denoise(Ypsd = psd_data$Ypsd, alpha = .01,
                        psd = psd, method = "fisherG")
Ypsd_den2 <- psd_denoise(Ypsd = psd_data$Ypsd, alpha = .01,
                         psd = psd, method = "BH")

bin_type <- "mean"
fbar <- binning(psd_data$freq, bin_size, bin_type = bin_type)
Ybar <- binning(psd_data$Ypsd, bin_size, bin_type = bin_type)
Ybar_den <- binning(Ypsd_den, bin_size, bin_type = bin_type)
Ybar_den2 <- binning(Ypsd_den2, bin_size, bin_type = bin_type)

par(mfrow = c(1,3))
xlim <- frange
type <- "l"
plot(fbar, Ybar, type = type, log = "xy",
     xlim = xlim, main = "No Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")
plot(fbar, Ybar_den, type = type, log = "xy",
     xlim = xlim, main = "fisherG Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")
plot(fbar, Ybar_den2, type = type, log = "xy",
     xlim = xlim, main = "BH Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")

# final fit
shof_obj$set_psd(freq = psd_data$freq, Ypsd = Ypsd_den)
shof_obj$set_est(est_type = "lp",
                 bin_size = bin_size, bin_type = "mean")

nlp <- shof_obj$nlp()
fit <- optim(par = fit$par, fn = nlp$fn, gr = nlp$gr, method = "BFGS")
lp_est <- to_est(fit$par, shof_obj)

# short version
set.seed(300)
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
lp_fit <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                  est_type = "lp", bin_size = bin_size)

# check
range(lp_fit$coef - lp_est$coef)
range(lp_fit$se - lp_est$se)
range((lp_fit$psd$Ypsd_fg - Ypsd_den)/Ypsd_den)

#' ## Fit MLE Estimator

# preliminary fit
set.seed(745)
# use LP for this
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
shof_obj$set_est(est_type = "lp",
                 bin_size = bin_size, bin_type = "median")
nlp <- shof_obj$nlp()
fit <- optim(par = phi0, fn = nlp$fn, gr = nlp$gr, method = "BFGS")

# denoising
theta_fit <- shof_obj$to_theta(phi = fit$par,
                               zeta = shof_obj$zeta(fit$par)$fn())
psd <- shof_psd(freq = psd_data$freq, theta = theta_fit, Temp = Temp)
Ypsd_den <- psd_denoise(Ypsd = psd_data$Ypsd, alpha = .01,
                        psd = psd, method = "fisherG")
Ypsd_den2 <- psd_denoise(Ypsd = psd_data$Ypsd, alpha = .01,
                         psd = psd, method = "BH")

bin_type <- "mean"
fbar <- binning(psd_data$freq, bin_size, bin_type = bin_type)
Ybar <- binning(psd_data$Ypsd, bin_size, bin_type = bin_type)
Ybar_den <- binning(Ypsd_den, bin_size, bin_type = bin_type)
Ybar_den2 <- binning(Ypsd_den2, bin_size, bin_type = bin_type)

par(mfrow = c(1,3))
xlim <- frange
type <- "l"
plot(fbar, Ybar, type = type, log = "xy",
     xlim = xlim, main = "No Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")
plot(fbar, Ybar_den, type = type, log = "xy",
     xlim = xlim, main = "fisherG Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")
plot(fbar, Ybar_den2, type = type, log = "xy",
     xlim = xlim, main = "BH Denoise")
lines(shof_obj$fbar,
      exp(shof_obj$zeta(fit$par)$fn()) * shof_obj$ufun(fit$par)$fn(),
      col = "red")

# final fit
shof_obj$set_psd(freq = psd_data$freq, Ypsd = Ypsd_den)
shof_obj$set_est(est_type = "mle",
                 bin_size = bin_size, bin_type = "mean")

nlp <- shof_obj$nlp()
fit <- optim(par = fit$par, fn = nlp$fn, gr = nlp$gr, method = "BFGS")
mle_est <- to_est(fit$par, shof_obj)

# short version
set.seed(745)
frange <- c(20, 60)*1000
psd_data <- psd_emp %>% filter(freq >= frange[1] & freq <= frange[2])
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
mle_fit <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                   est_type = "mle", bin_size = bin_size)

# check
range(mle_fit$coef - mle_est$coef)
range(mle_fit$se - mle_est$se)
range((mle_fit$psd$Ypsd_fg - Ypsd_den)/Ypsd_den)
