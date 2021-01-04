#' SHOF PSD with vector theta input.
#'
#' @param freq Frequency vector.
#' @param theta Parameter vector `theta = (k, f0, Q, Sf, alpha)`.
#' @param Temp Temperature
#' @return A vector of SHOF PSD values evaluated at the inputs.
shof_psd <- function(freq, theta, Temp) {
  theta <- as.list(theta)
  do.call(showf_psd,
          c(list(freq = freq, Temp = Temp, Sw = 0), theta))
  ## showf_psd(fseq, k, f0, Q, 0, Sf, alpha, Temp)
}

#' SHOW PSD with vector theta input.
#'
#' @param freq Frequency vector.
#' @param theta Parameter vector `theta = (k, f0, Q, Sw)`.
#' @param Temp Temperature
#' @return A vector of SHOW PSD values evaluated at the inputs.
show_psd <- function(freq, theta, Temp) {
  theta <- as.list(theta)
  do.call(showf_psd,
          c(list(freq = freq, Temp = Temp, Sf = 0, alpha = 0), theta))
  ## showf_psd(fseq, k, f0, Q, 0, Sf, alpha, Temp)
}

#' Convert a fitted value of `phi` to a complete estimate of `theta`.
#'
#' @param phi Vector of normalized PSD parameters.
#' @param obj PSD model object.
#' @return A list with elements `coef`, `vcov`, and `se`.
to_est <- function(phi, obj) {
  zeta <- obj$zeta()$fn(phi)
  coef <- obj$to_theta(phi, zeta)
  vcov <- obj$vcov(phi, zeta, to_theta = TRUE)
  colnames(vcov) <- rownames(vcov) <- names(coef)
  se <- sqrt(diag(vcov))
  list(coef = coef, se = se, vcov = vcov)
}

#' Fit PSD and return various statistics.
#'
#' @param psd_obj PSD model object.
#' @param theta0 Parameter values to initialize the optimization.
#' @param est_type Type of estimator: "nls", "lp", or "mle".
#' @param bin_size Bin size.
#' @param method Denoising method: "fisherG" or "BH".
#'
#' @return A list with elements:
#' - `psd`: A tibble with columns `freq`, `Ypsd`, `Ypsd_fg`, `Ypsd_bh`.
#' - `coef_pre`: Prefit parameter estimate on the `theta` scale.
#' - `coef`, `vcov`, `se`: All on the `theta` scale.
#'
#' @warning Modifies argument `psd_obj` externally (i.e. passes by reference).
fit_psd <- function(psd_obj, theta0,
                    est_type, bin_size = 100,
                    method = c("fisherG", "BH")) {
  method <- match.arg(method)
  # extract original data
  freq <- psd_obj$freq
  Ypsd <- psd_obj$Ypsd
  # preliminary fit
  psd_obj$set_est(est_type = ifelse(est_type == "mle", "lp", est_type),
                  bin_size = bin_size, bin_type = "median")
  phi0 <- psd_obj$to_phi(theta0)$phi
  nlp <- psd_obj$nlp()
  prefit <- optim(par = phi0, fn = nlp$fn, gr = nlp$gr, method = "BFGS")
  prefit$zeta <- psd_obj$zeta(prefit$par)$fn()
  # denoising
  theta_pre <- psd_obj$to_theta(phi = prefit$par, zeta = prefit$zeta)
  # to get psd of full frequency basis
  psd_obj$set_est(est_type = "mle", bin_type = "mean")
  psd_pre <- exp(prefit$zeta) * psd_obj$ufun(prefit$par)$fn()
  Ypsd_den <- psd_denoise(Ypsd = psd_obj$Ypsd, alpha = .01,
                          psd = psd_pre, method = method)
  ## Ypsd_fg <- psd_denoise(Ypsd = psd_obj$Ypsd, alpha = .01,
  ##                        psd = psd_pre, method = "fisherG")
  ## Ypsd_bh <- psd_denoise(Ypsd = psd_obj$Ypsd, alpha = .01,
  ##                        psd = psd_pre, method = "BH")
  # final fit
  psd_obj$set_psd(freq = freq, Ypsd = Ypsd_den)
  psd_obj$set_est(est_type = est_type,
                  bin_size = bin_size, bin_type = "mean")
  nlp <- psd_obj$nlp()
  fit <- optim(par = prefit$par, fn = nlp$fn, gr = nlp$gr, method = "BFGS")
  est <- to_est(fit$par, psd_obj)
  # reset psd in psd_obj
  psd_obj$set_psd(freq = freq, Ypsd = Ypsd)
  # output
  out <- list(psd = tibble(freq = freq,
                           Ypsd = Ypsd,
                           Ypsd_den = Ypsd_den),
              coef_pre = theta_pre)
  c(out, est)
}

#' Test plot for denoising.
#'
#' @param fit A list returned by [fit_psd()].
#' @param model The PSD model used: "shof" or "show".
#' @param Temp Temperature (Kelvin).
#' @param bin_size,bin_type How to bin for plotting.
plot_denoise <- function(fit, model = c("shof", "show"), Temp,
                         bin_size, bin_type) {
  # binning
  fbar <- binning(fit$psd$freq, bin_size = bin_size, bin_type = bin_type)
  Ybar <- binning(fit$psd$Ypsd, bin_size = bin_size, bin_type = bin_type)
  Ybar_den <- binning(fit$psd$Ypsd_den, bin_size, bin_type = bin_type)
  ## Ybar_fg <- binning(fit$psd$Ypsd_fg, bin_size, bin_type = bin_type)
  ## Ybar_bh <- binning(fit$psd$Ypsd_bh, bin_size, bin_type = bin_type)
  model <- match.arg(model)
  if(model == "shof") {
    psd_pre <- shof_psd(freq = fbar, theta = fit$coef_pre, Temp = Temp)
    psd_fit <- shof_psd(freq = fbar, theta = fit$coef, Temp = Temp)
  } else if(model == "show") {
    psd_pre <- show_psd(freq = fbar, theta = fit$coef_pre, Temp = Temp)
    psd_fit <- show_psd(freq = fbar, theta = fit$coef, Temp = Temp)
  }
  if(bin_type == "median") {
    Ybar <- Ybar/log(2)
    Ybar_den <- Ybar_den/log(2)
  }
  ## plot_data <- tibble(f = fbar, none = Ybar,
  ##                     fg = Ybar_fb, bh = Ybar_bh, fit = psd_fit)
  opar <- par()$mfrow
  par(mfrow = c(1,2))
  type <- "l"
  plot(fbar, Ybar, type = type, log = "xy", main = "No Denoise")
  lines(fbar, psd_pre, col = "red")
  legend("topright",
         legend = c("Original PSD", "Preliminary Fit"), fill = c("black", "red"))
  plot(fbar, Ybar_den, type = type, log = "xy", main = "fisherG Denoise")
  lines(fbar, psd_fit, col = "red")
  legend("topright",
         legend = c("Denoised PSD", "Final Fit"), fill = c("black", "red"))
  par(mfrow = opar)
}

