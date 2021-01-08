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

## ---- ou_make_psd_model, eval = FALSE-----------------------------------------
#  # create c++ file
#  ou_cpp <- make_psd_model(model = "OU",
#                           header = "OU_Model.hpp",
#                           class = "ou::UFun",
#                           ctor = "ou::make_Ufun",
#                           standalone = TRUE)

## ----ou_ufun, echo = FALSE, results = "asis"----------------------------------
cat("```cpp", readLines("OU_FitMethods.cpp"), "```", sep = "\n")

