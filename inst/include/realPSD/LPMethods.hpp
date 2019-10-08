/// @file LPMethods.hpp

#ifndef realPSD_LPMethods_hpp
#define realPSD_LPMethods_hpp 1

#include "LPFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type LP_zeta(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    // calculate zeta_hat
    return lp.zeta(Ubar);
  }

  template<class Type>
  Type LP_nlp(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    SIMULATE{
      // addtionally report the vector of residuals
      // FIXME: lp.zeta(Ubar) wastefully recomputes residuals
      matrix<Type> res(N,1);
      lp.res(res, Ubar);
      res.array() -= lp.zeta(Ubar);
      REPORT(res); 
    } 
    // calculate nlp
    return lp.nlp(Ubar);
  }

  template<class Type>
  Type LP_nll(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    PARAMETER(zeta);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    SIMULATE{
      // additionally report the vector of residuals
      matrix<Type> res(N,1);
      lp.res(res, Ubar);
      res.array() -= zeta;
      REPORT(res); 
    } 
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    Ubar = Ubar * fs;
    return lp.nll(Ubar, zeta);
  }

  
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
