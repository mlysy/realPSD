/// @file MLEMethods.hpp

#ifndef realPSD_MLEMethods_hpp
#define realPSD_MLEMethods_hpp 1

#include "MLEFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type MLE_methods(objective_function<Type>* obj,
		   const std::string& method,
		   const matrix<Type>& f,
		   const matrix<Type>& Y,
		   const matrix<Type>& U) {
    int N = f.size();
    MLE<Type> mle(N);
    mle.set_Y(Y);
    if(method == "MLE_tau") {
      return mle.tau(U);
    } else if(method == "MLE_nlp") {
      return mle.nlp(U);
    } else if(method == "MLE_nll") {
      PARAMETER(tau);
      return mle.nll(U, tau);
    } else {
      error("Unknown MLE method.");
    }
    return Type(0.0);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
