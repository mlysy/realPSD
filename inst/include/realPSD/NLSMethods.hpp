/// @file NLSMethods.hpp

#ifndef realPSD_NLSMethods_hpp
#define realPSD_NLSMethods_hpp 1

#include "NLSFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type NLS_methods(objective_function<Type>* obj,
		   const std::string& method,
		   const matrix<Type>& fbar,
		   const matrix<Type>& Ybar,
		   const matrix<Type>& Ubar) {
    int N = fbar.size();
    NLS<Type> nls(N);
    nls.set_Ybar(Ybar);
    if(method == "NLS_tau") {
      return nls.tau(Ubar);
    } else if(method == "NLS_nlp") {
      return nls.nlp(Ubar);
    } else if(method == "NLS_res") {
      // calculate residuals
      matrix<Type> res(N,1);
      nls.res(res, Ubar);
      ADREPORT(res); // make residuals differentiable
      return Type(0.0);
    } else if(method == "NLS_nll") {
      PARAMETER(tau);
      return nls.nll(Ubar, tau);
    } else {
      error("Unknown NLS method.");
    }
    return Type(0.0);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
