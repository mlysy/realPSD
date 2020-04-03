/// @file ComplexPoly.cpp

#include <Rcpp.h>
using namespace Rcpp;
#include <complex>

/// Evaluate a complex PSD-polynomial using only real operations.
///
/// @brief A PSD polynomial is defined as
///
/// ~~~
/// poly(x) = a_0 + a_1 * (1i * x) + a_2 * (1i * x)^2 + ... + a_p * (1i *x)^p,
/// ~~~
///
/// where `a_0, ..., a_p` and `x` are real.
void complex_poly(double& y_real, double& y_imag,
		  const double x, const double* alpha, const int degree) {
  int n_odd = (degree + 1) / 2; // number of odd powers
  bool last_even = (degree % 2) == 0; // whether there is one more even power
  int sgn = 2 * (degree % 4 < 2) - 1; // sign of largest even power
  double x2 = x*x;
  // initialize
  y_real = sgn * alpha[degree - !last_even];
  y_imag = n_odd > 0 ? sgn * alpha[degree - last_even] : 0.0;
  for(int ii=degree-2; ii>0; ii-=2) {
    sgn = -sgn;
    y_real = sgn * alpha[ii - !last_even] + x2 * y_real;
    y_imag = sgn * alpha[ii - last_even] + x2 * y_imag;
  }
  // last even power
  if(last_even && degree > 0) y_real = alpha[0] + x2 * y_real;
  y_imag *= sgn * x; // last_even and last_odd have different leading signs
  return;
}

// [[Rcpp::export]]
std::complex<double> cpoly(double x, NumericVector alpha) {
  double y_real, y_imag;
  complex_poly(y_real, y_imag, x, REAL(alpha), alpha.length()-1);
  std::complex<double> y(y_real, y_imag);
  return y;
}
