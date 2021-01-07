## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE, 
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
# link to packages
pkg_link <- function(pkg, link) {
  if(link == "github") {
    link <- paste0("https://github.com/mlysy/", pkg)
  } else if(link == "cran") {
    link <- paste0("https://CRAN.R-project.org/package=", pkg)
  }
  paste0("[**", pkg, "**](", link, ")")
}
cran_link <- function(pkg) pkg_link(pkg, "cran")
github_link <- function(pkg) pkg_link(pkg, "github")

## ---- SHOW_R6, eval = FALSE---------------------------------------------------
#  model <- sample(c("SHOW_nat", "SHOW_log", "SHOW_comp"), size = 1) # can be chosen from any of these
#  your_model <- R6::R6Class(
#            classname = "your_model",
#            inherit = realPSD::psd_model,
#            private = list(
#              Temp_ = NULL, # temperature
#              n_phi_ = 3, # specify the number of parameters, for SHOW model, f0, Q, Rw
#              tmb_DLL_ = "realPSD_TMBExports",
#              tmb_model_ = model
#            ),
#            active = list(
#              #' @field Temp Temperature of the system.
#              Temp = function() private$Temp_
#            )
#          )

## ---- SHOW_MLE_ufun, eval = FALSE---------------------------------------------
#  # get MLE Ufun
#  scale <- mean(Ypsd)
#  # create TMB model object
#  tmod <- TMB::MakeADFun(data = list(model = model,
#                                     method = "MLE_ufun",
#                                     f = matrix(freq),
#                                     Y = matrix(Ypsd/scale),
#                                     scale = log(scale)),
#                         parameters = list(phi = matrix(rep(0, 3))),
#                         ADreport = TRUE,
#                         silent = TRUE, DLL = "realPSD_TMBExports")

## ---- ou_make_psd_model, eval = FALSE-----------------------------------------
#  # create c++ file
#  ou_cpp <- make_psd_model(model = "OU",
#                           header = "OU_Model.hpp",
#                           class = "ou::UFun",
#                           ctor = "ou::make_Ufun",
#                           standalone = TRUE)

## ----ou_ufun, echo = FALSE, results = "asis"----------------------------------
cat("```cpp", readLines("OU_FitMethods.cpp"), "```", sep = "\n")

