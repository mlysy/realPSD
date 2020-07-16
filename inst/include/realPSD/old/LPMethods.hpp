/// @file LPMethods.hpp

#ifndef realPSD_LPMethods_hpp
#define realPSD_LPMethods_hpp 1

#include "LPFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type LP_methods(objective_function<Type>* obj,
		  const std::string& method,
		  const matrix<Type>& fbar,
		  const matrix<Type>& Zbar,
		  const matrix<Type>& Ubar) {
    int N = fbar.size();
    LP<Type> lp(N);
    lp.set_Zbar(Zbar);
    if(method == "LP_zeta") {
      return lp.zeta(Ubar);
    } else if(method == "LP_nlp") {
      return lp.nlp(Ubar);
    } else if(method == "LP_res") {
      matrix<Type> res(N,1); // output variable
      lp.res(res, Ubar); // calculate residuals
      ADREPORT(res); // set them to be differentiable
      return Type(0.0);
    } else if(method == "LP_nll") {
      PARAMETER(zeta);
      return lp.nll(Ubar, zeta);
    } else {
      // std::string msg = "Unknown method " + method + ".";
      error("Unknown LP method.");
    }
    return Type(0.0);
  }		  
  
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
