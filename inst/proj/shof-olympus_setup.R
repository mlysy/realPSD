## ---- shof_setup
## source(system.file("proj", "plot_psd.R", package = "realPSD"))
## if(use_local) {
##   source("shof-functions.R")
## } else {
##   source(system.file("proj", "shof-functions.R", package = "realPSD"))
##   ## source(system.file("tests", "testthat", "realPSD-testfunctions.R",
##   ##                    package = "realPSD"))
## }
source(file.path(proj_path, "shof-functions.R"))
## tmb_recompile()
## require(R6)

#' Pure R implementation of the SHOF normalized PSD.
shof_ufun_r <- function(freq, phi) {
  ephi <- exp(phi)
  f2 <- (freq/ephi[1])^2
  ephi[3]/freq^phi[4] + 1/((1-f2)^2 + f2/ephi[2]^2)
}

#' Pure R implementation of the SHOF normalized PSD.
show_ufun_r <- function(freq, phi) {
  ephi <- exp(phi)
  f2 <- (freq/ephi[1])^2
  ephi[3] + 1/((1-f2)^2 + f2/ephi[2]^2)
}

#' Standard errors plot.
#'
#' @param par_labels Named vector of strings.  Each is the display expression of a parameter, named as the factor level in `theta`.
#' @param est_labels Named vector of strings; same as `par_labels` for estimator names (`nls`, `lp`, `mle`).
#' @param col Color for estimators.
#' @param lgd_pos Position of legend (ntc units relative to top figure).
#' @param scale_breaks `breaks` argument for y-axis.
#' @param ylab_size Size of y-axis labels.
se_plot <- function(sho_fit, par_labels, est_labels, col, lgd_pos,
                    scale_breaks, ylab_size) {
  ## plt_est <- shof_fit[!is.na(shof_fit)] %>%
  line_size <- .6
  plt_est <- sho_fit %>%
    mutate(theta = factor(theta, levels = names(par_labels),
                          labels = par_labels),
           est_type = factor(est_type, level = names(est_labels),
                             labels = est_labels)) %>%
    ggplot(aes(x = bin_size, y = coef, group = est_type)) +
    geom_line(aes(color = est_type), size = line_size) +
    scale_color_manual(values = col, guide = FALSE) +
    scale_y_continuous(breaks = scale_breaks) +
    ## scale_y_log10() +
    facet_wrap(~ theta, scales = "free", labeller = "label_parsed") +
    ylab("(a) Estimate") +
    xlab(NULL) +
    labs(col = "Estimator") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = ylab_size),
      strip.background = element_blank(),
      strip.placement = "outside",
      ## legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.justification = "left",
      ## legend.box.background = element_blank()
      ## legend.box.background = element_rect(size = .4)
      ## legend.position = lgd_pos,
      ## legend.justification = c("right", "top"),
      legend.box.background = element_rect(size = .2),
      ## legend.box.just = "right"
    )
  plt_se <- sho_fit %>%
    mutate(theta = factor(theta, levels = names(par_labels),
                          labels = par_labels),
           est_type = factor(est_type, level = names(est_labels),
                             labels = est_labels)) %>%
    ggplot(aes(x = bin_size, y = se, group = est_type)) +
    geom_line(aes(color = est_type), size = line_size) +
    geom_hline(yintercept = 0, alpha = 0.25, linetype = "dotted") +
    scale_color_manual(values = col, guide = FALSE) +
    ## scale_y_log10() +
    scale_y_continuous(breaks = scale_breaks) +
    facet_wrap(~ theta, scales = "free", labeller = "label_parsed") +
    ylab("(b) Standard Error") +
    xlab(NULL) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = ylab_size),
      strip.background = element_blank(),
      ## strip.text.x = element_blank(),
      strip.placement = "outside")
  plt_cv <- sho_fit %>%
    mutate(theta = factor(theta, levels = names(par_labels),
                          labels = par_labels),
           est_type = factor(est_type, level = names(est_labels),
                             labels = est_labels)) %>%
    ggplot(aes(x = bin_size, y = se/coef, group = est_type)) +
    geom_line(aes(color = est_type), size = line_size) +
    geom_hline(yintercept = 0, alpha = 0.25, linetype = "dotted") +
    ## geom_hline(yintercept = 0) +
    scale_color_manual(values = col) +
    ## scale_y_log10() +
    scale_y_continuous(breaks = scale_breaks) +
    facet_wrap(~ theta, scales = "free", labeller = "label_parsed") +
    ## ylab("(c) Coefficient of Variation") +
    ylab("(c) CV") +
    xlab("\nBin Size") +
    labs(col = "Estimator") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = ylab_size),
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.background = element_rect(size = .2)
    )
  list(est = plt_est, se = plt_se, cv = plt_cv)
}
