#' ---
#' title: "Application: Calibration of an Atomic Force Microscope"
#' ---

#' ## Prerequisites

#+ shof_setup
# packages
require(realPSD)
require(tidyverse)
require(optimCheck)
require(numDeriv)
require(gridExtra)
require(RColorBrewer)

# extra source code
# use current dir for proj folder.  otherwise installed package's.
proj_local <- TRUE
save_data <- FALSE # save data
if(proj_local) {
  proj_path <- getwd()
} else {
  proj_path <- system.file("proj", package = "realPSD")
}
source(file.path(proj_path, "shof-olympus_setup.R"))

#' ## Plot PSD Data

#+ data
calc_psd <- FALSE
fs <- 5e6 # sampling frequency (Hz)
Temp <- 298 # temperature (Kelvin)
frng_eig1 <- c(20, 60)*1000
if(!calc_psd) {
  ## stop("haven't done this yet")
  ## psd_emp <- readRDS("olympus-psd_sine.rds")
  psd_bin <- readRDS(file = "olympus-psd_bin.rds")
  psd_data <- readRDS(file = "olympus-psd_data.rds")
} else {
  Xt <- readRDS("olympus-time_series_noisy.rds")
  ## scale_factor <- 1
  psd_emp <- as_tibble(periodogram(Xt, fs = fs))
  ## %>% mutate(Ypsd = Ypsd/(2*pi * freq * 20)^2)
  psd_emp$Ypsd <- psd_emp$Ypsd/(2*pi * psd_emp$freq * 20)^2
  # binning (for plotting)
  psd_bin <- sapply(psd_emp, binning,
                    bin_size = 100, bin_type = "mean") %>%
    as_tibble()
  # psd data for fitting
  psd_data <- psd_emp %>% filter(freq >= frng_eig1[1] & freq <= frng_eig1[2])
  if(save_data) {
    saveRDS(psd_bin, file = "olympus-psd_bin.rds")
    saveRDS(psd_data, file = "olympus-psd_data.rds")
  }
}

# plot data
## clrs <- c(brewer.pal(3, "Set1")[2], "black")
## psd_color <- "deepskyblue1"
psd_color <- brewer.pal(3, "Set1")[1]
tick_breaks <- 10^(-30:30)
minor_breaks <- as.matrix(expand.grid(1:9, tick_breaks))
minor_breaks <- minor_breaks[,1] * minor_breaks[,2]
label10 <- function(x) parse(text = paste0("10^", log10(x)))
line_size <- .3
## psd_bin <- readRDS("psd_bin.rds")
plt_main <- psd_bin %>%
  transmute(x = freq/1000, y = Ypsd) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(size = line_size, color = psd_color) +
  scale_x_log10(
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
  ## ylab(expression("PSD "*(m^2/Hz))) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("(a) Periodogram") +
  annotation_logticks(size = 0.3) +
  theme_bw()
## plt_main

frange <- frng_eig1
plt_eig1 <- psd_bin %>%
  filter(frange[1] <= freq & freq <= frange[2]) %>%
  transmute(x = freq/1000, y = Ypsd, mode = 1)
plt_eig1 <- plt_eig1 %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(size = line_size, color = psd_color) +
  ## scale_x_log10() +
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
## plt_eig1

frange <- c(100, 300)*1000
plt_eig2 <- psd_bin %>%
  filter(frange[1] <= freq & freq <= frange[2]) %>%
  transmute(x = freq/1000, y = Ypsd)
plt_eig2 <- plt_eig2 %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(size = line_size, color = psd_color) +
  ## scale_x_log10() +
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
## plt_eig2

plot_group <- grid.arrange(
  grobs = list(plt_main, plt_eig1, plt_eig2),
  layout_matrix = rbind(c(1,1), c(2,3)),
  left = grid::textGrob(expression("PSD "*(m^2/Hz)),
                        rot = 90),
  bottom = grid::textGrob(expression("Frequency (KHz)"))
)
plot_group

#' ## Fit the SHOF Model

#' ### Test `shof_model`

#+ test_psd
# plausible parameters
theta0 <- c(k = .178, f0 = 33563.03, Q = 52.8,
            Sf = 7.218499e-18,
            alpha = 1.99)

# compare to TMB
bin_size <- 100
shof_obj <- shof_model$new(freq = psd_data$freq, Ypsd = psd_data$Ypsd,
                           Temp = Temp, bin_size = bin_size)

phi0 <- shof_obj$to_phi(theta0)$phi
zeta0 <- shof_obj$to_phi(theta0)$zeta
fbar <- binning(psd_data$freq, bin_size = bin_size, bin_type = "mean")
U_tmb <- shof_obj$ufun()$fn(phi0)
U_r <- shof_ufun_r(freq = fbar, phi = phi0)
range(U_tmb - U_r)

## plot(fbar, U_tmb, type = "l", log = "xy")
## lines(fbar, U_r, type = "l", col = "red")


#' ### Fit the NLS Estimator

shof_fit <- setNames(vector("list", 3), c("nls", "lp", "mle"))

#+ nls_fit
## frange <- frng_eig1
## psd_data <- psd_emp %>% filter(freq >= frange[1] & freq <= frange[2])
shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
shof_fit$nls <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                        est_type = "nls", bin_size = bin_size)

# plot denoise
plot_denoise(fit = shof_fit$nls, Temp = Temp, model = "shof",
             bin_size = bin_size, bin_type = "mean")

## saveRDS(nls_fit, file = "shof-nls_fit.rds")

#' ## Fit LP estimator

#+ lp_fit
## shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
shof_fit$lp <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                       est_type = "lp", bin_size = bin_size)

# plot denoise
plot_denoise(fit = shof_fit$lp, Temp = Temp, model = "shof",
             bin_size = bin_size, bin_type = "mean")
## plot_denoise(psd = shof_fit$lp$psd,
##              theta_pre = shof_lp_fit$coef_pre, theta_fit = shof_lp_fit$coef,
##              bin_size = bin_size, bin_type = "mean")

## saveRDS(shof_lp_fit, file = "shof-lp_fit.rds")

#' ## Fit MLE estimator

#+ mle_fit
## shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
shof_fit$mle <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                        est_type = "mle", bin_size = bin_size)

# plot denoise
plot_denoise(fit = shof_fit$mle, Temp = Temp, model = "shof",
             bin_size = bin_size, bin_type = "mean")
## plot_denoise(psd = mle_fit$psd,
##              theta_pre = mle_fit$coef_pre, theta_fit = mle_fit$coef,
##              bin_size = bin_size, bin_type = "mean")

## saveRDS(mle_fit, file = "shof-mle_fit.rds")

if(save_data) {
  saveRDS(shof_fit, file = "olympus-shof_fit.rds")
}

#' ## Fit the SHOW estimators

#+ show_test
# initial values
theta0 <- c(k = .178, f0 = 33563.03, Q = 52.8, Sw = 5e-27)
show_obj <- show_model$new(freq = psd_data$freq, Ypsd = psd_data$Ypsd,
                           Temp = Temp)

# compare TMB to R
bin_size <- 100
show_obj$set_est(bin_size = bin_size)
phi0 <- show_obj$to_phi(theta0)$phi
zeta0 <- show_obj$to_phi(theta0)$zeta
fbar <- binning(psd_data$freq, bin_size = bin_size, bin_type = "mean")
U_tmb <- show_obj$ufun()$fn(phi0)
U_r <- show_ufun_r(freq = fbar, phi = phi0)
range(U_tmb - U_r)

plot(show_obj$fbar, show_obj$Ybar, type = "l", log = "xy")
lines(show_obj$fbar, show_psd(show_obj$fbar, theta0, Temp), col = "red")

#+ show_fit
show_fit <- setNames(vector("list", 3), c("nls", "lp", "mle"))

# nls
show_fit$nls <- fit_psd(psd_obj = show_obj, theta0 = theta0,
                        est_type = "nls", bin_size = bin_size)

# plot denoise
plot_denoise(fit = show_fit$nls, Temp = Temp, model = "show",
             bin_size = bin_size, bin_type = "mean")


# lp
show_fit$lp <- fit_psd(psd_obj = show_obj, theta0 = theta0,
                       est_type = "lp", bin_size = bin_size)

plot_denoise(fit = show_fit$lp, Temp = Temp, model = "show",
             bin_size = bin_size, bin_type = "mean")


# mle (with lp prefit)
# FIXME: this fails
if(FALSE) {
  show_fit$mle <- fit_psd(psd_obj = show_obj, theta0 = theta0,
                          est_type = "mle", bin_size = bin_size)
}

if(save_data) {
  saveRDS(show_fit, file = "olympus-show_fit.rds")
}


#' ## Compare Estimates
#'
#' ### Tabular Form

#+ fit_tab
# look at the fits
tnames <- c("f0", "k", "Q")
sho_coef <- rbind(shof_nls = shof_fit$nls$coef[tnames],
                  shof_lp = shof_fit$lp$coef[tnames],
                  shof_mle = shof_fit$mle$coef[tnames],
                  show_nls = show_fit$nls$coef[tnames],
                  show_lp = show_fit$lp$coef[tnames])
sho_coef[,1] <- sho_coef[,1]/1000
sho_se <- rbind(shof_nls = shof_fit$nls$se[tnames],
                shof_lp = shof_fit$lp$se[tnames],
                shof_mle = shof_fit$mle$se[tnames],
                show_nls = show_fit$nls$se[tnames],
                show_lp = show_fit$lp$se[tnames])
sho_se[,1] <- sho_se[,1]/1000
(abs(shof_fit$nls$coef - shof_fit$lp$coef)/shof_fit$lp$se)[tnames]

digits <- c(5, 3, 3)
sapply(c(f0 = 1, k = 2, Q = 3), function(ii) {
  paste0(signif(sho_coef[,ii], digits[ii]), " (",
         signif(sho_se[,ii], 2), ")")
})

signif(sho_coef, 4)
signif(sho_se, 2)

#' ### Graphical Form

#+ fit_plot
bin_size <- 100
plot_data <- tibble(
  fbar = binning(psd_data$freq, bin_size = bin_size, bin_type = "mean"),
  Ybar = binning(psd_data$Ypsd, bin_size = bin_size, bin_type = "mean"),
  Ybar_den = binning(shof_fit$lp$psd$Ypsd_den, bin_size = bin_size,
                     bin_type = "median")/log(2),
  shof_nls = shof_psd(freq = fbar, theta = shof_fit$nls$coef, Temp = Temp),
  shof_lp = shof_psd(freq = fbar, theta = shof_fit$lp$coef, Temp = Temp),
  shof_mle = shof_psd(freq = fbar, theta = shof_fit$mle$coef, Temp = Temp),
  show_nls = show_psd(freq = fbar, theta = show_fit$nls$coef, Temp = Temp),
  show_lp = show_psd(freq = fbar, theta = show_fit$lp$coef, Temp = Temp)
)

tick_breaks <- 10^(-30:30)
minor_breaks <- as.matrix(expand.grid(1:9, tick_breaks))
minor_breaks <- minor_breaks[,1] * minor_breaks[,2]
label10 <- function(x) parse(text = paste0("10^", log10(x)))
labels <- c("PSD - Original", "PSD - Denoised",
            "SHOW - NLS", "SHOW - LP",
            "SHOF - NLS", "SHOF - LP", "SHOF - MLE", " ", "  ")
levels <- c("Ybar", "Ybar_den",
            "show_nls", "show_lp",
            "shof_nls", "shof_lp", "shof_mle", "blank1", "blank2")
## lgd_order <- c(1:2, 8, 5:7, 3:4, 9)
lgd_order <- 1:8
lgd_nrow <- 4
## lgd_order <- 1:9
## names(labels) <- levels
labels <- labels[1:length(lgd_order)]
levels <- levels[1:length(lgd_order)]
breaks <- labels[lgd_order] # labels[c(1:2, 5:7, 3:4)]
nfit <- length(labels) - 2
clrs <- brewer.pal(4, "Set1")
clrs <- c("grey80",
          "black",
          ## "black",
          ## clrs[c(4, 4, 1:3)],
          rep("orange", 2),
          ## clrs[1:2],
          clrs[1:3],
          rep("white", 2))
clrs <- clrs[lgd_order]
linetype <- c(rep("solid", 3),
              rep("dashed", 1),
              rep("solid", 4),
              rep("solid", 2))
## linetype <- c(rep("solid", 2), rep("dashed", 2), rep("solid", 3))
linetype <- linetype[lgd_order]
line_size <- c(.4, .4,
               .7, .7,
               ## 1, 1,
               1, 1, 1,
               0, 0)
line_size <- line_size[lgd_order]
drop <- FALSE
plot_data %>%
  pivot_longer(cols = Ybar:show_lp,
               names_to = "type", values_to = "Ybar") %>%
  mutate(type = factor(type,
                       levels = levels,
                       labels = labels)) %>%
ggplot(aes(x = fbar/1000, y = Ybar, group = type)) +
  geom_line(aes(color = type, size = type, linetype = type)) +
  scale_color_manual(values = clrs, breaks = breaks, drop = drop) +
  scale_size_manual(values = line_size, breaks = breaks, drop = drop) +
  scale_linetype_manual(values = linetype, breaks = breaks, drop = drop) +
  scale_y_log10(minor_breaks = minor_breaks,
                labels = label10) +
  annotation_logticks(size = 0.3, sides = "l") +
  xlab(expression("Frequency (kHz)")) +
  ylab(expression("PSD "*(m^2/Hz))) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.box.background = element_rect(size = .3),
        legend.box.just = "right") +
  guides(color = guide_legend(nrow = lgd_nrow),
         size = guide_legend(nrow = lgd_nrow),
         linetype = guide_legend(nrow = lgd_nrow))

#' ## SHOF with Varying Bin Size

#+ shof_bin
theta0 <- c(k = .178, f0 = 33563.03, Q = 52.8,
            Sf = 7.218499e-18, alpha = 1.99)
job_descr <- expand.grid(est_type = c("nls", "lp", "mle"),
                         bin_size = seq(10, 1000, by = 10),
                         method = "fisherG",
                         stringsAsFactors = FALSE)
## job_descr <- rbind(job_descr, c(est_type = "mle", bin_size = 100, method = "BH"))
njobs <- nrow(job_descr)
calc_bin <- FALSE
if(!calc_bin) {
  shof_bin <- readRDS("olympus-shof_bin.rds")
} else {
  shof_obj <- shof_model$new(freq = psd_data$freq, Ypsd = psd_data$Ypsd,
                             Temp = Temp)
  ## shof_obj2 <- shof_model$new(freq = psd_data$freq, Ypsd = psd_data$Ypsd,
  ##                             Temp = Temp)
  ## identical(shof_obj$Ypsd, shof_obj2$Ypsd)
  # use same mle fit for all sets
  shof_fit_mle <- fit_psd(psd_obj = shof_obj, theta0 = theta0,
                          est_type = "mle", bin_size = 100)
  job_fun <- function(ijob) {
    message("ijob = ", ijob)
    bin_size <- job_descr$bin_size[ijob]
    est_type <- job_descr$est_type[ijob]
    method <- job_descr$method[ijob]
    if(est_type == "mle") {
      fit <- shof_fit_mle
    } else {
      fit <- tryCatch(
        fit_psd(psd_obj = shof_obj, theta0 = theta0,
                est_type = est_type, bin_size = bin_size, method = method),
        error = function(e) NULL
      )
      if(is.null(fit)) {
        # reset psd
        shof_obj$set_psd(freq = psd_data$freq, Ypsd = psd_data$Ypsd)
        return(NA)
      }
    }
    tibble(est_type = est_type, bin_size = bin_size,
           theta = names(fit$coef),
           coef = as.numeric(fit$coef),
           se = as.numeric(fit$se))
  }
  # nls and lp fit
  shof_fit <- lapply(1:njobs, job_fun)
  ## saveRDS(shof_fit, file = "olympus-shof_bin.rds")
}

# summary for plotting
sho_fit <- shof_bin[!is.na(shof_bin)] %>%
  bind_rows() %>%
  filter(theta %in% c("f0", "k", "Q")) %>%
  # convert to kHz
  pivot_wider(names_from = theta, values_from = coef:se) %>%
  mutate(coef_f0 = coef_f0 / 1000,
         se_f0 = se_f0 / 1000) %>%
  pivot_longer(cols = coef_k:se_Q,
               names_to = c("stat", "theta"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = stat, values_from = value)


source(file.path(proj_path, "shof-olympus_setup.R"))
par_labels <- c(f0 = "f[0]*'  '*(kHz)",
                k = "k*'  '*(N/m)", Q = "Q*'  (unitless)'")
est_labels <- c(nls = "NLS", lp = "LP", mle = "MLE")
est_clrs <- brewer.pal(3, "Set1")
lgd_pos <- c(.14, .80)
## scale_breaks <- waiver()
scale_breaks <- pretty
ylab_size <- 10
plt <- se_plot(sho_fit = sho_fit,
               par_labels = par_labels, est_labels = est_labels,
               col = est_clrs, lgd_pos = lgd_pos,
               scale_breaks = scale_breaks, ylab_size = ylab_size)
grid::grid.newpage()
grid::grid.draw(rbind(ggplotGrob(plt$est),
                      ggplotGrob(plt$se),
                      ggplotGrob(plt$cv)))

## sho_fit %>%
##   group_by(theta) %>%
##   summarize(coef_min = min(coef),
##             coef_max = max(coef),
##             se_min = min(se),
##             se_max = max(se))

## sho_fit %>%
##   filter(theta == "f0") %>%
##   ggplot(aes(x = bin_size, y = coef, group = est_type)) +
##   geom_line(aes(color = est_type))
