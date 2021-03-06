# test generic fou model header created by make_psd_model
require(realPSD)
require(whisker)
require(TMB)
require(TMBtools)
require(Rcpp)
require(testthat)

# standalone == FALSE, the output file should be put under a package built by TMBtools
make_psd_model(
  name = "fou",
  header = "fou.hpp",
  class = "fou::UFun",
  ctor = "fou::make_Ufun",
  include = "realPSD/FitMethods.hpp",
  method = "realPSD::FitMethods",
  standalone = FALSE,
  path = "fou_Generics.hpp"
)

tmb_create_package(path = "../../../testPackage",
    tmb_files = "fou_Generics.hpp")
getwd()

# Then we need to put "fou.hpp" to the same directory of fou_Generics.hpp
# and change the DESCRIPTION file of the testPackage created by TMBtools linking to realPSD
# (to make sure the header file "realPSD/FitMethods.hpp" can be found)
# We can successfully compile the build testPackage

