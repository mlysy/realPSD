/// @file NLSMethods.hpp

#ifndef realPSD_NLSMethods_hpp
#define realPSD_NLSMethods_hpp 1

#include "NLSFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type NLS_tau(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Ybar);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    NLS<Type> nls(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    // calculate tau_hat
    return nls.tau(Ubar);
  }

  template<class Type>
  Type NLS_nll(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Ybar);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    PARAMETER(tau);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    NLS<Type> nls(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    // FIXME: shouldn't wastefully compute nll if only residuals are desired.
    SIMULATE{
      // calculate the vector of residuals
      matrix<Type> res(N,1);
      // res = nls.res(Ubar);
      nls.res(res, Ubar, tau);
      REPORT(res); 
    } 
    // calculate nll
    return nls.nll(Ubar, tau);
  }

  template<class Type>
  Type NLS_nlp(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Ybar);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    NLS<Type> nls(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    // FIXME: calculate only one of nlp or res.
    SIMULATE{
      // calculate the vector of residuals
      matrix<Type> res(N,1);
      Type tau = nls.tau(Ubar);
      // res = nls.res(Ubar);
      nls.res(res, Ubar, tau);
      REPORT(res); 
    } 
    // calculate nlp
    return nls.nlp(Ubar);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
