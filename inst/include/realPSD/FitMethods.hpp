/// @file FitMethods.hpp

#ifndef realPSD_FitMethods_hpp
#define realPSD_FitMethods_hpp 1

#include "LPMethods.hpp"
#include "NLSMethods.hpp"
#include "MLEMethods.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type, class UFun>
  Type FitMethods(objective_function<Type>* obj) {
    // pick method
    DATA_STRING(method);
    if(method == "UFun") {
      // data
      DATA_MATRIX(f);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate U
      int N = f.size();
      UFun Ufun(N);
      Ufun.set_f(f);
      SIMULATE {
	matrix<Type> U(N,1);
	Ufun.eval(U, phi);
	REPORT(U);
      }
      return Type(0);
    } else if(method.find("LP_") == 0) {
      // data
      DATA_MATRIX(fbar);
      DATA_MATRIX(Zbar);
      DATA_VECTOR(fs);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate Ubar
      int N = fbar.size();
      matrix<Type> Ubar(N,1);
      UFun Ufun(N);
      Ufun.set_f(fbar);
      Ufun.eval(Ubar, phi);
      Ubar = Ubar * fs;
      return LP_methods(obj, method, fbar, Zbar, Ubar);
    } else if(method.find("MLE_") == 0) {
      // data
      DATA_MATRIX(f);
      DATA_MATRIX(Y);
      DATA_VECTOR(fs);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate U
      int N = f.size();
      matrix<Type> U(N,1);
      UFun Ufun(N);
      Ufun.set_f(f);
      Ufun.eval(U, phi);
      U = U * fs;
      return MLE_methods(obj, method, f, Y, U);
    } else if(method.find("NLS_") == 0) {
      // data
      DATA_MATRIX(fbar);
      DATA_MATRIX(Ybar);
      DATA_VECTOR(fs);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate Ubar
      int N = fbar.size();
      matrix<Type> Ubar(N,1);
      UFun Ufun(N);
      Ufun.set_f(fbar);
      Ufun.eval(Ubar, phi);
      Ubar = Ubar * fs;
      return NLS_methods(obj, method, fbar, Ybar, Ubar);
      error("Unknown method.");
    }
    return Type(0.0);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
