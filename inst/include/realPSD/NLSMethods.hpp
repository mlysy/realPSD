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
    // FIXME: add residual method
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
    // addtionally report the vector of residuals
    // calculate the residuals
    // FIXME: RES should be res
    SIMULATE{
      matrix<Type> RES(N,1);
      RES = nls.res(Ubar);
      REPORT(RES); 
    } 
    // calculate nlp
    return nls.nlp(Ubar);
  }

  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
