/// @file ComplexPoly.cpp

#include <Rcpp.h>
using namespace Rcpp;

/// Evaluate a complex polynomial using only real operations.
void complex_poly(double& y_real, double& y_imag,
		  const double x, const double* alpha, const int degree) {
  int n = degree + 1;
  int n_odd = n / 2; // number of odd terms
  bool extra
  y_real = 0.0;
  y_imag = 0.0;
}
