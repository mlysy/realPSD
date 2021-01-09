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
# load required packages
require(realPSD)
require(R6)
# get current working directory
wd <- getwd()
# invisibly copy OU_Model.hpp from realPSD to the current working directory
# this trick helps realPSD pass win_builder test in rebuilding vignettes section
file.copy(from = system.file("include", "realPSD", "OU_Model.hpp", 
                              package = "realPSD"),
          to = wd)
# create c++ file
ou_cpp <- make_psd_model(model = "OU",
                         header = "OU_Model.hpp",
                         class = "ou::UFun",
                         ctor = "ou::make_Ufun",
                         standalone = TRUE)
# write the content into a .cpp file called OU_FitMethods.cpp under the current working directory
cat(ou_cpp, sep = "\n", file = "OU_FitMethods.cpp")

## ---- ou_make_psd_model, eval = FALSE-----------------------------------------
#  # create c++ file
#  ou_cpp <- make_psd_model(model = "OU",
#                           header = "OU_Model.hpp",
#                           class = "ou::UFun",
#                           ctor = "ou::make_Ufun",
#                           standalone = TRUE)

## ----ou_ufun, echo = FALSE, results = "asis"----------------------------------
cat("```cpp", readLines("OU_FitMethods.cpp"), "```", sep = "\n")

