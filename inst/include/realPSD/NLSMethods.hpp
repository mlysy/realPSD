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
    // calculate Ubar
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
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
    // calculate Ubar
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
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
    // calculate Ubar
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    // calculate nlp
    return nls.nlp(Ubar);
  }

  template<class Type>
  Type NLS_res(objective_function<Type>* obj) {
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
    // output variables
    matrix<Type> res(N,1);
    // calculate Ubar
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    // calculate residuals
    nls.res(res, Ubar);
    ADREPORT(res);
    return Type(0.0);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
